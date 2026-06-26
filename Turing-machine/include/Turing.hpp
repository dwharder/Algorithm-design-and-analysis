#ifndef CA_UWATERLOO_DWHARDER_TURING_MACHINE_LIBRARY
#define CA_UWATERLOO_DWHARDER_TURING_MACHINE_LIBRARY

// This implementation was developed through discussion with,
// and with significant assistance from, ChatGPT.

#include <cstddef>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>
#include <concepts>
#include <iostream>
#include <sstream>
#include <utility>
#include <regex>
#include <iomanip>

// Class and struct declarations
struct State;

template <typename State_type = std::size_t>
class Run_result;

template <typename Symbol>
class Alphabet;

template <typename Symbol>
class Tape;

template <typename Symbol>
class Turing_program;

// Operator declarations
std::ostream &operator<<( std::ostream &out, State s );

template <typename State_type>
std::ostream &operator<<( std::ostream &out, Run_result<State_type> const &r );

template <typename Symbol>
std::ostream &operator<<( std::ostream &out, Tape<Symbol> const &tape );

// ================================================================
// TM constants
// ================================================================

namespace TM {
    using state_type = std::size_t;
    using code_type  = std::size_t;

    inline constexpr state_type HALT =
        std::numeric_limits<state_type>::max();

    inline constexpr state_type TRUE  = HALT - 1;
    inline constexpr state_type FALSE = HALT - 2;
    inline constexpr state_type FAIL  = HALT - 3;

    struct Blank_t {};
    struct All_t {};
    struct Same_t {};

    inline constexpr Blank_t BLANK{};
    inline constexpr All_t   ALL{};
    inline constexpr Same_t  SAME{};

    enum class Move { left, right, stay };

    inline constexpr Move LEFT  = Move::left;
    inline constexpr Move RIGHT = Move::right;
    inline constexpr Move STAY  = Move::stay;
}

inline std::string state_to_string( TM::state_type s ) {
    if ( s == TM::HALT )  return "HALT";
    if ( s == TM::TRUE )  return "TRUE";
    if ( s == TM::FALSE ) return "FALSE";
    if ( s == TM::FAIL )  return "FAIL";

    return std::to_string( s );
}

template <typename Symbol>
using initial_symbol_type =
    std::variant<Symbol, TM::Blank_t>;

  //////////////////
 // struct State //
//////////////////

struct State {
    TM::state_type value;
};

inline std::ostream &operator<<( std::ostream &out, State s ) {
    return out << state_to_string( s.value );
}

  //////////////////////////////////
 // class Run_result<State_type> //
//////////////////////////////////

template <typename State_type>
class Run_result {
public:
    explicit constexpr Run_result( State_type state )
        : state_{ state } {}

    constexpr operator State_type() const {
        return state_;
    }

    State_type state() const {
        return state_;
    }

private:
    State_type state_;
};

template <typename State_type>
std::ostream &operator<<( std::ostream &out, Run_result<State_type> const &r ) {
    return out << State{ r.state() };
}

  ////////////////////////////
 // class Alphabet<Symbol> //
////////////////////////////

template <typename Symbol>
class Alphabet {
public:
    using symbol_type = Symbol;
    using code_type   = TM::code_type;

    Alphabet( std::initializer_list<Symbol> symbols ) {
        for ( Symbol const &s : symbols ) {
            insert( s );
        }
    }

    Alphabet( std::string_view symbols )
    requires std::same_as<Symbol, char> {
        for ( char c : symbols ) {
            insert( c );
        }
    }

    explicit Alphabet( code_type n )
    requires std::is_integral_v<Symbol> {
        for ( code_type k{ 0 }; k < n; ++k ) {
            insert( static_cast<Symbol>( k ) );
        }
    }

    Alphabet( Symbol first, Symbol last )
    requires std::is_integral_v<Symbol> {
        if ( first > last ) {
            throw std::invalid_argument{ "Invalid alphabet range." };
        }

        // It may happen that 'last' is such that ++last overflows.
        for ( Symbol s{ first }; true; ++s ) {
            insert( s );

            if ( s == last ) {
                break;
            }
        }
    }

    code_type encode( Symbol const &s ) const {
        auto it{ forward_.find( s ) };

        if ( it == forward_.end() ) {
            throw std::out_of_range{ "Symbol not in alphabet." };
        }

        return it->second;
    }

    Symbol decode( code_type code ) const {
        if ( code >= reverse_.size() ) {
            throw std::out_of_range{ "Symbol code outside alphabet." };
        }

        return reverse_[code];
    }

    code_type size() const {
        return reverse_.size();
    }

    static constexpr code_type BLANK =
        std::numeric_limits<code_type>::max();

    static constexpr bool is_blank( code_type code ) {
        return code == BLANK;
    }

    std::string code_to_string( code_type code ) const {
        if ( is_blank( code ) ) {
            return "BLANK";
        }

        Symbol symbol{ decode( code ) };

        if constexpr ( std::same_as<Symbol, char> ) {
            return std::string{ 1, symbol };
        } else if constexpr ( std::is_arithmetic_v<Symbol> ) {
            return std::to_string( symbol );
        } else {
            std::ostringstream out;
            out << symbol;
            return out.str();
        }
    }

    template <typename Function>
    void for_each_symbol( Function f ) const {
        for ( Symbol const &s : reverse_ ) {
            f( s );
        }
    }

private:
    std::unordered_map<Symbol, code_type> forward_;
    std::vector<Symbol> reverse_;

    void insert( Symbol const &s ) {
        if constexpr ( std::same_as<Symbol, char> ) {
            if ( s == '_' ) {
                throw std::invalid_argument(
                    "The '_' character is reserved for displaying blank cells."
                );
            }
        }

        if ( forward_.contains( s ) ) {
            throw std::invalid_argument{
               "Duplicate symbol are not allowed within an alphabet"
            };
        }

        code_type code{ reverse_.size() };
        forward_[s] = code;
        reverse_.push_back( s );
    }
};

// Class template argument deduction guides
template <typename Symbol>
Alphabet( std::initializer_list<Symbol> ) -> Alphabet<Symbol>;

template <typename Symbol>
Alphabet( Symbol ) -> Alphabet<Symbol>;

template <typename Symbol>
Alphabet( Symbol, Symbol ) -> Alphabet<Symbol>;

Alphabet( std::string_view ) -> Alphabet<char>;
Alphabet( char const * ) -> Alphabet<char>;

  ////////////////////////
 // class Tape<Symbol> //
////////////////////////
template <typename Symbol>
class Tape {
public:
    using alphabet_type  = Alphabet<Symbol>;
    using symbol_type    = Symbol;
    using code_type      = std::size_t;
    using position_type  = std::ptrdiff_t;

    static constexpr code_type BLANK =
        std::numeric_limits<code_type>::max();

    explicit Tape(
        alphabet_type const &alphabet,
        position_type initial_head = 0
    ) :
        alphabet_{ alphabet },
        head_{ initial_head },
        initial_head_{ initial_head },
        min_touched_{ initial_head },
        max_touched_{ initial_head }
    {
        // Intentionally empty.
    }

    Tape(
        alphabet_type const &alphabet,
        std::initializer_list<initial_symbol_type<Symbol>> symbols,
        position_type initial_head = 0
    ) :
        alphabet_{ alphabet },
        head_{ initial_head },
        initial_head_{ initial_head },
        min_touched_{ initial_head },
        max_touched_{ initial_head }
    {
        position_type pos{ 0 };

        for ( auto const &entry : symbols ) {
            if ( std::holds_alternative<Symbol>( entry ) ) {
                write_at( pos, std::get<Symbol>( entry ) );
            } else {
                note_access( pos );   // important for embedded blanks
            }

            ++pos;
        }
    }

    Tape(
        alphabet_type const &alphabet,
        std::string_view symbols,
        position_type initial_head = 0
    ) requires std::same_as<Symbol, char> :
    Tape{ alphabet, initial_head }
    {
        position_type pos{ 0 };

        for ( char c : symbols ) {
            write_at( pos, c );
            ++pos;
        }
    }

    Tape(
        alphabet_type const &alphabet,
        std::initializer_list<std::string_view> words,
        position_type initial_head = 0
    ) requires std::same_as<Symbol, char>
        : alphabet_{ alphabet },
          head_{ initial_head },
          initial_head_{ initial_head },
          min_touched_{ initial_head },
          max_touched_{ initial_head }
    {
        position_type pos{ 0 };
        bool first_word{ true };

        for ( std::string_view word : words ) {
            if ( !first_word ) {
                note_access( pos );  // the separating blank
                ++pos;
            }

            for ( char c : word ) {
                write_at( pos, c );
                ++pos;
            }

            first_word = false;
        }
    }

    alphabet_type const &alphabet() const {
        return alphabet_;
    }

    position_type head() const {
        return head_;
    }

    code_type read_code() const {
        return read_code_at( head_ );
    }

    Symbol read() const {
        code_type code{ read_code() };

        if ( code == BLANK ) {
            throw std::logic_error{
                "Cannot decode a blank symbol."
            };
        }

        return alphabet_.decode( code );
    }

    bool is_blank() const {
        return read_code() == BLANK;
    }

    void write( Symbol const &symbol ) {
        write_code_at( head_, alphabet_.encode( symbol ) );
    }

    void write_code( code_type code ) {
        validate_code( code );
        write_code_at( head_, code );
    }

    void erase() {
        write_code_at( head_, BLANK );
    }

    void move_left() {
        --head_;
        note_access( head_ );
    }

    void move_right() {
        ++head_;
        note_access( head_ );
    }

    void stay() {
        note_access( head_ );
    }

    void reset_head() {
        head_ = initial_head_;
        note_access( head_ );
    }

    code_type read_code_at( position_type pos ) const {
        auto it{ cells_.find( pos ) };

        if ( it == cells_.end() ) {
            return BLANK;
        }

        return it->second;
    }

    void write_at( position_type pos, Symbol const &symbol ) {
        write_code_at( pos, alphabet_.encode( symbol ) );
    }

    void erase_at( position_type pos ) {
        write_code_at( pos, BLANK );
    }

    position_type min_touched() const {
        return min_touched_;
    }

    position_type max_touched() const {
        return max_touched_;
    }

    void print(
        std::ostream &out = std::cout,
        position_type padding = 3,
        bool show_head = true
    ) const {
        position_type first{ min_touched_ - padding };
        position_type last{ max_touched_ + padding };

        for ( position_type pos{ first }; pos <= last; ++pos ) {
            code_type code{ read_code_at( pos ) };

            if ( show_head && (pos == head()) ) out << "[";

            if ( code == BLANK ) {
                out << '_';
            } else {
                out << alphabet().code_to_string( code );
            }

            if ( show_head && (pos == head()) ) out << "]";

            if ( pos != last ) {
                out << ' ';
            }
        }

        out << '\n';
    }

private:
    // Member variables //
    alphabet_type const &alphabet_;

    std::unordered_map<position_type, code_type> cells_;

    position_type head_;
    position_type initial_head_;

    position_type min_touched_;
    position_type max_touched_;

    // Private helper functions //
    void validate_code( code_type code ) const {
        if ( code == BLANK ) {
            return;
        }

        if ( code >= alphabet_.size() ) {
            throw std::out_of_range{
                "The symbol code is outside the alphabet."
            };
        }
    }

    void write_code_at( position_type pos, code_type code ) {
        validate_code( code );
        note_access( pos );

        if ( code == BLANK ) {
            cells_.erase( pos );
        } else {
            cells_[pos] = code;
        }
    }

    void note_access( position_type pos ) {
        if ( pos < min_touched_ ) {
            min_touched_ = pos;
        }

        if ( pos > max_touched_ ) {
            max_touched_ = pos;
        }
    }
};

// Class template argument deduction guides

template <typename Symbol>
Tape( Alphabet<Symbol> const &, std::ptrdiff_t )
    -> Tape<Symbol>;

template <typename Symbol>
Tape(
    Alphabet<Symbol> const &,
    std::initializer_list<Symbol>,
    std::ptrdiff_t
) -> Tape<Symbol>;

Tape(
    Alphabet<char> const &,
    std::string_view,
    std::ptrdiff_t
) -> Tape<char>;

Tape(
    Alphabet<char> const &,
    char const *,
    std::ptrdiff_t
) -> Tape<char>;

Tape(
    Alphabet<char> const &,
    std::initializer_list<std::string_view>
) -> Tape<char>;

Tape(
    Alphabet<char> const &,
    std::initializer_list<std::string_view>,
    std::ptrdiff_t
) -> Tape<char>;

// Operator definitions

template <typename Symbol>
std::ostream &operator<<( std::ostream &out, Tape<Symbol> const &tape ) {
    using position_type = typename Tape<Symbol>::position_type;
    using code_type     = typename Tape<Symbol>::code_type;

    out << "~~ ";

    // Print tape
    for ( position_type pos = tape.min_touched(); pos <= tape.max_touched(); ++pos ) {
        code_type code = tape.read_code_at( pos );

        if ( pos == tape.head() ) out << "[";

        if ( code == Tape<Symbol>::BLANK ) {
            out << "_";
        } else {
            out << tape.alphabet().code_to_string( code );
        }

        if ( pos == tape.head() ) out << "]";

        out << " ";
    }

    return out << "~~";
}

  //////////////////////////////////
 // class Turing_program<Symbol> //
//////////////////////////////////

template <typename Symbol>
class Turing_program {
public:
    using alphabet_type = Alphabet<Symbol>;
    using tape_type     = Tape<Symbol>;

    using state_type = TM::state_type;
    using code_type  = TM::code_type;
    using Move       = TM::Move;

    static constexpr state_type INITIAL_STATE = 0;

private:
    struct Key {
        state_type state;
        code_type  symbol;

        bool operator==( Key const &other ) const {
            return state == other.state && symbol == other.symbol;
        }
    };

    struct Key_hash {
        std::size_t operator()( Key const &key ) const {
            return std::hash<state_type>{}( key.state )
                 ^ (std::hash<code_type>{}( key.symbol ) << 1 );
        }
    };

public:
    struct Rule_rhs {
        std::variant<state_type, TM::Same_t> new_state;
        std::variant<Symbol, TM::Blank_t, TM::Same_t> new_symbol;
        Move move;

        Rule_rhs( state_type s, Symbol const &sym, Move m )
            : new_state{ s }, new_symbol{ sym }, move{ m } {}

        Rule_rhs( state_type s, TM::Blank_t b, Move m )
            : new_state{ s }, new_symbol{ b }, move{ m } {}

        Rule_rhs( state_type s, TM::Same_t same, Move m )
            : new_state{ s }, new_symbol{ same }, move{ m } {}

        Rule_rhs( TM::Same_t same_state, Symbol const &sym, Move m )
            : new_state{ same_state }, new_symbol{ sym }, move{ m } {}

        Rule_rhs( TM::Same_t same_state, TM::Blank_t b, Move m )
            : new_state{ same_state }, new_symbol{ b }, move{ m } {}

        Rule_rhs( TM::Same_t same_state, TM::Same_t same_symbol, Move m )
            : new_state{ same_state }, new_symbol{ same_symbol }, move{ m } {}

        Rule_rhs(
            std::variant<state_type, TM::Same_t> s,
            std::variant<Symbol, TM::Blank_t, TM::Same_t> sym,
            Move m
        )
            : new_state{ std::move( s ) },
              new_symbol{ std::move( sym ) },
              move{ m }
        {}
    };

private:
    struct Step_transition {
        state_type new_state;
        code_type  new_symbol;
        Move       move;
    };

public:
    struct Call_transition {
        Turing_program const *program;
        std::unordered_map<state_type, state_type> return_states;
    };

private:
    using Transition =
        std::variant<Step_transition, Call_transition>;

    using State_selector =
        std::vector<state_type>;

    using Symbol_selector =
    std::variant<
        std::vector<code_type>,
        TM::All_t,
        alphabet_type const *
    >;

public:
    class Transition_proxy {
    public:
        Transition_proxy(
            Turing_program &program,
            State_selector states,
            Symbol_selector symbols,
            bool allow_override
        )
            : program_{ program },
              states_{ std::move( states ) },
              symbols_{ std::move( symbols ) },
              allow_override_{ allow_override }
        {
            for ( state_type state : states_ ) {
                program_.validate_input_state( state );
            }
        }

        Transition_proxy &operator=( Rule_rhs const &rhs ) {
            if ( std::holds_alternative<state_type>( rhs.new_state ) ) {
                program_.validate_output_state(
                    std::get<state_type>( rhs.new_state )
                );
            }

            if (
                std::holds_alternative<TM::All_t>( symbols_ )
                && std::holds_alternative<TM::Same_t>( rhs.new_symbol )
            ) {
                for ( state_type s : states_ ) {
                    bool same_state{
                        std::holds_alternative<TM::Same_t>( rhs.new_state )
                        || std::get<state_type>( rhs.new_state ) == s
                    };

                    if ( same_state ) {
                        std::clog
                            << "Warning: state "
                            << state_to_string( s )
                            << " has a TM::ALL transition that preserves both "
                               "the state and the symbol. Since TM::ALL includes "
                               "blanks, this transition has no stopping condition.\n";
                    }
                }
            }

            for_each_expanded_symbol(
                [&]( state_type in_state, code_type in_symbol ) {
                    state_type out_state{
                        std::holds_alternative<TM::Same_t>( rhs.new_state )
                            ? in_state
                            : std::get<state_type>( rhs.new_state )
                    };

                    code_type out_symbol{
                        std::holds_alternative<TM::Same_t>( rhs.new_symbol )
                            ? in_symbol
                            : program_.encode_output_symbol( rhs.new_symbol )
                    };

                    program_.insert_transition(
                        Key{ in_state, in_symbol },
                        Step_transition{
                            out_state,
                            out_symbol,
                            rhs.move
                        },
                        allow_override_
                    );
                }
            );

            return *this;
        }

        Transition_proxy &operator=( Call_transition const &call ) {
            if ( call.program == nullptr ) {
                throw std::invalid_argument{
                    "A call transition cannot call a null program."
                };
            }

            for_each_expanded_symbol(
                [&]( state_type in_state, code_type in_symbol ) {
                    program_.insert_transition(
                        Key{ in_state, in_symbol },
                        call,
                        allow_override_
                    );
                }
            );

            return *this;
        }

    private:
        Turing_program &program_;
        State_selector states_;
        Symbol_selector symbols_;
        bool allow_override_;

        template <typename Function>
        void for_each_expanded_symbol( Function f ) const {
            for ( state_type s : states_ ) {
                if ( std::holds_alternative<TM::All_t>( symbols_ ) ) {
                    for ( code_type k{ 0 }; k < program_.alphabet_.size(); ++k ) {
                        f( s, k );
                    }

                    f( s, tape_type::BLANK );

                } else if ( std::holds_alternative<alphabet_type const *>( symbols_ ) ) {
                    auto const *alphabet{ std::get<alphabet_type const *>( symbols_ ) };

                    if ( alphabet != &program_.alphabet_ ) {
                        throw std::invalid_argument{
                            "The alphabet selector must be the program's alphabet."
                        };
                    }

                    for ( code_type k{ 0 }; k < program_.alphabet_.size(); ++k ) {
                        f( s, k );
                    }

                } else {
                    auto const &codes{
                        std::get<std::vector<code_type>>( symbols_ )
                    };

                    for ( code_type c : codes ) {
                        f( s, c );
                    }
                }
            }
        }
    };

    explicit Turing_program(
        alphabet_type const &alphabet,
        std::initializer_list<state_type> terminal_states = { TM::HALT }
    )
        : alphabet_{ alphabet }
    {
        if ( terminal_states.size() == 0 ) {
            throw std::invalid_argument{
                "A Turing program must have at least one halting state."
            };
        }

        for ( state_type s : terminal_states ) {
            terminal_states_.insert( s );
        }
    }

    // Character alphabets

    Turing_program(
        alphabet_type const &alphabet,
        std::initializer_list<std::string_view> rules,
        std::initializer_list<state_type> halting_states = { TM::HALT }
    ) requires std::same_as<Symbol, char>
        : Turing_program{ alphabet, halting_states }
    {
        validate_rule_alphabet_for_char();
        parse_rules( rules );
    }

    Turing_program(
        std::string name,
        alphabet_type const &alphabet,
        std::initializer_list<std::string_view> rules,
        std::initializer_list<state_type> halting_states = { TM::HALT }
    ) requires std::same_as<Symbol, char>
        : Turing_program{ alphabet, rules, halting_states }
    {
        register_program( std::move( name ) );
    }

    // Integer alphabets

    Turing_program(
        alphabet_type const &alphabet,
        std::initializer_list<std::string_view> rules,
        std::initializer_list<state_type> halting_states = { TM::HALT }
    ) requires ( std::is_integral_v<Symbol> && !std::same_as<Symbol, char> )
        : Turing_program{ alphabet, halting_states }
    {
        validate_rule_alphabet_for_integer();
        parse_rules( rules );
    }

    Turing_program(
        std::string name,
        alphabet_type const &alphabet,
        std::initializer_list<std::string_view> rules,
        std::initializer_list<state_type> halting_states = { TM::HALT }
    ) requires ( std::is_integral_v<Symbol> && !std::same_as<Symbol, char> )
        : Turing_program{ alphabet, rules, halting_states }
    {
        register_program( std::move( name ) );
    }

    // ----- ordinary rule definition -----

    Transition_proxy operator()( state_type state, Symbol const &symbol ) {
        return make_proxy(
            { state },
            symbols_to_codes( { symbol } ),
            false
        );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        Symbol const &symbol
    ) {
        return make_proxy(
            State_selector{ states },
            symbols_to_codes( { symbol } ),
            false
        );
    }

    Transition_proxy operator()(
        state_type state,
        std::initializer_list<Symbol> symbols
    ) {
        return make_proxy(
            { state },
            symbols_to_codes( symbols ),
            false
        );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        std::initializer_list<Symbol> symbols
    ) {
        return make_proxy(
            State_selector{ states },
            symbols_to_codes( symbols ),
            false
        );
    }

    Transition_proxy operator()( state_type state, TM::Blank_t ) {
        return make_proxy(
            { state },
            std::vector<code_type>{ tape_type::BLANK },
            false
        );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        TM::Blank_t
    ) {
        return make_proxy(
            State_selector{ states },
            std::vector<code_type>{ tape_type::BLANK },
            false
        );
    }

    Transition_proxy operator()( state_type state, TM::All_t ) {
        return make_proxy( { state }, TM::ALL, false );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        TM::All_t
    ) {
        return make_proxy( State_selector{ states }, TM::ALL, false );
    }

    Transition_proxy operator()( state_type state, alphabet_type const &alphabet ) {
        return make_proxy( { state }, &alphabet, false );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        alphabet_type const &alphabet
    ) {
        return make_proxy( State_selector{ states }, &alphabet, false );
    }

    Transition_proxy operator()( state_type state, std::string_view symbols )
    requires std::same_as<Symbol, char> {
        return make_proxy(
            { state },
            symbols_to_codes( symbols ),
            false
        );
    }

    Transition_proxy operator()(
        std::initializer_list<state_type> states,
        std::string_view symbols
    ) requires std::same_as<Symbol, char> {
        return make_proxy(
            State_selector{ states },
            symbols_to_codes( symbols ),
            false
        );
    }

    // ----- explicit override rule definition -----

    Transition_proxy override_rule( state_type state, Symbol const &symbol ) {
        return make_proxy(
            { state },
            symbols_to_codes( { symbol } ),
            true
        );
    }

    Transition_proxy override_rule( state_type state, TM::Blank_t ) {
        return make_proxy(
            { state },
            std::vector<code_type>{ tape_type::BLANK },
            true
        );
    }

    Transition_proxy override_rule( state_type state, TM::All_t ) {
        return make_proxy( { state }, TM::ALL, true );
    }

    Transition_proxy override_rule(
        state_type state,
        std::initializer_list<Symbol> symbols
    ) {
        return make_proxy(
            { state },
            symbols_to_codes( symbols ),
            true
        );
    }

    Transition_proxy override_rule(
        std::initializer_list<state_type> states,
        std::initializer_list<Symbol> symbols
    ) {
        return make_proxy(
            State_selector{ states },
            symbols_to_codes( symbols ),
            true
        );
    }

    Transition_proxy override_rule(
        state_type state,
        alphabet_type const &alphabet
    ) {
        return make_proxy( { state }, &alphabet, true );
    }

    Transition_proxy override_rule(
        std::initializer_list<state_type> states,
        alphabet_type const &alphabet
    ) {
        return make_proxy( State_selector{ states }, &alphabet, true );
    }

    Transition_proxy override_rule(
        state_type state,
        std::string_view symbols
    ) requires std::same_as<Symbol, char> {
        return make_proxy(
            { state },
            symbols_to_codes( symbols ),
            true
        );
    }

    Transition_proxy override_rule(
        std::initializer_list<state_type> states,
        std::string_view symbols
    ) requires std::same_as<Symbol, char> {
        return make_proxy(
            State_selector{ states },
            symbols_to_codes( symbols ),
            true
        );
    }

    // ----- calls -----

    Call_transition call(
        Turing_program const &program,
        std::initializer_list<
            std::pair<state_type const, state_type>
        > returns = {}
    ) const {
        return Call_transition{
            &program,
            std::unordered_map<state_type, state_type>{ returns }
        };
    }

    // ----- execution -----

    Run_result<state_type> operator()(
        tape_type &tape,
        state_type start_state = INITIAL_STATE,
        std::size_t max_steps = 100000
    ) const {
        return run( tape, start_state, max_steps );
    }

    Run_result<state_type> run(
        tape_type &tape,
        state_type start_state = INITIAL_STATE,
        std::size_t max_steps = 100000
    ) const {
        if ( &tape.alphabet() != &alphabet_ ) {
            throw std::invalid_argument{
                "The tape and program do not use the same alphabet."
            };
        }

        state_type state{ start_state };

        for ( std::size_t k{ 0 }; k < max_steps; ++k ) {
            if ( is_terminal_state( state ) ) {
                return Run_result<state_type>{ state };
            }

            state = execute_one_transition( tape, state, max_steps );
        }

        throw std::runtime_error{
            "The program did not terminate within the maximum number of steps."
        };
    }

    void add_terminal_state( state_type state ) {
        terminal_states_.insert( state );
    }

    bool is_terminal_state( state_type state ) const {
        return terminal_states_.contains( state );
    }

    alphabet_type const &alphabet() const {
        return alphabet_;
    }

    Run_result<state_type> trace(
        tape_type &tape,
        std::ostream &out = std::cout,
        state_type start_state = INITIAL_STATE,
        std::size_t max_steps = 100000
    ) const {
        if ( &tape.alphabet() != &alphabet_ ) {
            throw std::invalid_argument{
                "The tape and program do not use the same alphabet."
            };
        }

        tape_type dry_run_tape{ tape };
        run( dry_run_tape, start_state, max_steps );

        auto first{ dry_run_tape.min_touched() };
        auto last{ dry_run_tape.max_touched() };
        auto width{ state_width() };

        state_type state{ start_state };

        for ( std::size_t k{ 0 }; k < max_steps; ++k ) {
            out << std::setw( width )
                << state_to_string( state )
                << "  ";

            print_trace_tape( out, tape, first, last );
            out << '\n';

            if ( is_terminal_state( state ) ) {
                return Run_result<state_type>{ state };
            }

            state = execute_one_transition( tape, state, max_steps );
        }

        throw std::runtime_error{
            "The program did not terminate within the maximum number of steps."
        };
    }

private:
    alphabet_type const &alphabet_;
    std::unordered_map<Key, Transition, Key_hash> transitions_;
    std::unordered_set<state_type> terminal_states_;

    std::string name_;

    inline static std::unordered_map<
        std::string,
        Turing_program const *
    > registry_;

    Transition_proxy make_proxy(
        State_selector states,
        Symbol_selector symbols,
        bool allow_override
    ) {
        return Transition_proxy{
            *this,
            std::move( states ),
            std::move( symbols ),
            allow_override
        };
    }

    std::vector<code_type> symbols_to_codes(
        std::initializer_list<Symbol> symbols
    ) const {
        std::vector<code_type> result;

        for ( Symbol const &s : symbols ) {
            result.push_back( alphabet_.encode( s ) );
        }

        return result;
    }

    std::vector<code_type> symbols_to_codes( std::string_view symbols ) const
    requires std::same_as<Symbol, char> {
        std::vector<code_type> result;

        for ( char c : symbols ) {
            result.push_back( alphabet_.encode( c ) );
        }

        return result;
    }

    code_type encode_output_symbol(
        std::variant<Symbol, TM::Blank_t, TM::Same_t> const &symbol
    ) const {
        if ( std::holds_alternative<TM::Blank_t>( symbol ) ) {
            return tape_type::BLANK;
        }

        if ( std::holds_alternative<TM::Same_t>( symbol ) ) {
            throw std::logic_error{
                "TM::SAME cannot be encoded without an input symbol."
            };
        }

        return alphabet_.encode( std::get<Symbol>( symbol ) );
    }

    void insert_transition(
        Key key,
        Transition transition,
        bool allow_override
    ) {
        if ( allow_override ) {
            transitions_[key] = std::move( transition );
            return;
        }

        auto [it, inserted]{
            transitions_.emplace( key, std::move( transition ) )
        };

        if ( !inserted ) {
            throw std::logic_error{
                "A transition is already defined for this state and symbol."
            };
        }
    }

    state_type execute_one_transition(
        tape_type &tape,

        state_type state,
        std::size_t max_steps
    ) const {
        code_type symbol{ tape.read_code() };

        auto it{ transitions_.find( Key{ state, symbol } ) };

        if ( it == transitions_.end() ) {
            std::ostringstream message;

            message << "No transition is defined for state "
                    << state_to_string( state )
                    << " and symbol "
                    << alphabet_.code_to_string( symbol )
                    << ".";

            throw std::logic_error{ message.str() };
        }

        Transition const &transition{ it->second };

        if ( std::holds_alternative<Step_transition>( transition ) ) {
            Step_transition const &t{
                std::get<Step_transition>( transition )
            };

            tape.write_code( t.new_symbol );

            switch ( t.move ) {
            case Move::left:
                tape.move_left();
                break;

            case Move::right:
                tape.move_right();
                break;

            case Move::stay:
                tape.stay();
                break;
            }

            return t.new_state;
        }

        Call_transition const &call{
            std::get<Call_transition>( transition )
        };

        state_type result{
            static_cast<state_type>(
                call.program->run( tape, INITIAL_STATE, max_steps )
            )
        };

        auto mapped{ call.return_states.find( result ) };

        if ( mapped != call.return_states.end() ) {
            return mapped->second;
        }

        if ( is_terminal_state( result ) ) {
            return result;
        }

        throw std::runtime_error{
            "The called program returned an unmapped halting state "
            + state_to_string( result )
            + " that is not a halting state of the calling program."
        };
    }

    static bool is_predefined_halting_state( state_type state ) {
        return state == TM::HALT
            || state == TM::TRUE
            || state == TM::FALSE
            || state == TM::FAIL;
    }

    void validate_input_state( state_type state ) const {
        if ( is_terminal_state( state ) ) {
            throw std::logic_error{
                "A transition cannot be defined from a halting state."
            };
        }

        if ( is_predefined_halting_state( state ) ) {
            throw std::logic_error{
                "The pre-defined halting state " + state_to_string( state ) + " was not declared in this program."
            };
        }
    }

    void validate_output_state( state_type state ) const {
        if (
            is_predefined_halting_state( state )
            && !is_terminal_state( state )
        ) {
            throw std::logic_error{
                "A transition cannot go to a predefined state unless "
                "that state was declared as a halting state."
            };
        }
    }

    void validate_rule_alphabet_for_char() const
    requires std::same_as<Symbol, char> {
        alphabet_.for_each_symbol(
            []( char c ) {
                if ( c == '*' || c == '#' || c == '=' ) {
                    throw std::invalid_argument{
                        "The characters '*', '#', and '=' are reserved "
                        "in string-based transition rules."
                    };
                }
            }
        );
    }

    void validate_rule_alphabet_for_integer() const
    requires ( std::is_integral_v<Symbol> && !std::same_as<Symbol, char> ) {
        alphabet_.for_each_symbol(
            []( Symbol s ) {
                if constexpr ( std::is_signed_v<Symbol> ) {
                    if ( s < 0 ) {
                        throw std::invalid_argument{
                            "String-based transition rules require non-negative "
                            "integer alphabet symbols."
                        };
                    }
                }
            }
        );
    }

    static state_type parse_state_token( std::string const &token ) {
        if ( token == "Z" ) return TM::HALT;
        if ( token == "T" ) return TM::TRUE;
        if ( token == "F" ) return TM::FALSE;
        if ( token == "X" ) return TM::FAIL;

        std::size_t used{};
        unsigned long long value{ std::stoull( token, &used ) };

        if ( used != token.size() ) {
            throw std::invalid_argument{ "Invalid state token: " + token };
        }

        return static_cast<state_type>( value );
    }

    static Move parse_move_token( std::string const &token ) {
        if ( token == "L" ) return TM::LEFT;
        if ( token == "R" ) return TM::RIGHT;
        if ( token == "S" ) return TM::STAY;

        throw std::invalid_argument{ "Invalid move token: " + token };
    }

    std::variant<state_type, TM::Same_t>
    parse_output_state_token( std::string const &token ) const {
        if ( token == "=" ) {
            return TM::SAME;
        }

        return parse_state_token( token );
    }

    State_selector parse_state_selector( std::string const &token ) const {
        if ( token.size() >= 2 && token.front() == '{' && token.back() == '}' ) {
            State_selector states;

            std::string body{ token.substr( 1, token.size() - 2 ) };
            std::istringstream in{ body };
            std::string part;

            while ( std::getline( in, part, ',' ) ) {
                if ( part.empty() ) {
                    throw std::invalid_argument{
                        "Empty state in state set: " + token
                    };
                }

                states.push_back( parse_state_token( part ) );
            }

            return states;
        }

        return State_selector{ parse_state_token( token ) };
    }

    Symbol_selector parse_symbol_selector( std::string const &token ) const {
        if ( token == "*" ) {
            return TM::ALL;
        }

        if ( token == "#" ) {
            return &alphabet_;
        }

        if ( token == "_" ) {
            return std::vector<code_type>{ tape_type::BLANK };
        }

        if ( token.size() >= 2 && token.front() == '{' && token.back() == '}' ) {
            std::vector<code_type> codes;

            std::string body{ token.substr( 1, token.size() - 2 ) };
            std::istringstream in{ body };
            std::string part;

            while ( std::getline( in, part, ',' ) ) {
                if ( part.empty() ) {
                    throw std::invalid_argument{
                        "Empty symbol in symbol set: " + token
                    };
                }

                if ( part == "_" ) {
                    codes.push_back( tape_type::BLANK );
                } else {
                    codes.push_back(
                        alphabet_.encode( parse_symbol_token( part ) )
                    );
                }
            }

            return codes;
        }

        return std::vector<code_type>{
            alphabet_.encode( parse_symbol_token( token ) )
        };
    }

    std::variant<Symbol, TM::Blank_t, TM::Same_t>
    parse_output_symbol_token( std::string const &token ) const {
        if ( token == "=" ) {
            return TM::SAME;
        }

        if ( token == "_" ) {
            return TM::BLANK;
        }

        return parse_symbol_token( token );
    }

    Symbol parse_symbol_token( std::string const &token ) const
    requires ( std::is_integral_v<Symbol> && !std::same_as<Symbol, char> ) {
        std::size_t used{};
        long long value{ std::stoll( token, &used ) };

        if ( used != token.size() || value < 0 ) {
            throw std::invalid_argument{ "Invalid integer symbol: " + token };
        }

        if constexpr ( std::is_signed_v<Symbol> ) {
            if ( value > std::numeric_limits<Symbol>::max() ) {
                throw std::out_of_range{ "Integer symbol too large: " + token };
            }
        } else {
            using unsigned_long_long = unsigned long long;

            auto unsigned_value{ static_cast<unsigned_long_long>( value ) };

            if ( unsigned_value > std::numeric_limits<Symbol>::max() ) {
                throw std::out_of_range{ "Integer symbol too large: " + token };
            }
        }

        return static_cast<Symbol>( value );
    }

    Symbol parse_symbol_token( std::string const &token ) const
    requires std::same_as<Symbol, char> {
        if ( token.size() != 1 ) {
            throw std::invalid_argument{
                "Invalid character symbol: " + token
            };
        }

        return token[0];
    }

    std::pair<state_type const, state_type>
    parse_return_mapping( std::string const &token ) const {
        auto pos{ token.find( "->" ) };

        if ( pos == std::string::npos ) {
            throw std::invalid_argument{
                "Invalid return-state mapping: " + token
            };
        }

        std::string from{ token.substr( 0, pos ) };
        std::string to{ token.substr( pos + 2 ) };

        if ( from.empty() || to.empty() ) {
            throw std::invalid_argument{
                "Invalid return-state mapping: " + token
            };
        }

        return {
            parse_state_token( from ),
            parse_state_token( to )
        };
    }

    void parse_rules( std::initializer_list<std::string_view> rules ) {
        for ( std::string_view rule : rules ) {
            parse_rule( rule );
        }
    }

    void parse_rule( std::string_view rule ) {
        std::istringstream in{ std::string{ rule } };

        std::string in_state_token;
        std::string in_symbol_token;
        std::string arrow;
        std::string out_state_token;

        if ( !( in >> in_state_token >> in_symbol_token >> arrow >> out_state_token ) ) {
            throw std::invalid_argument{
                "Invalid transition rule: " + std::string{ rule }
            };
        }

        if ( arrow != "->" ) {
            throw std::invalid_argument{
                "Expected '->' in transition rule: " + std::string{ rule }
            };
        }

        State_selector in_states{ parse_state_selector( in_state_token ) };

        if ( out_state_token == "call" ) {
            std::string program_name;

            if ( !( in >> program_name ) ) {
                throw std::invalid_argument{
                    "Expected program name after 'call' in rule: " + std::string{ rule }
                };
            }

            Turing_program const &called_program{ find_program( program_name ) };

            std::vector<std::pair<state_type const, state_type>> mappings;
            std::string mapping_token;

            while ( in >> mapping_token ) {
                mappings.push_back( parse_return_mapping( mapping_token ) );
            }

            Call_transition call_transition{
                &called_program,
                std::unordered_map<state_type, state_type>{
                    mappings.begin(),
                    mappings.end()
                }
            };

            assign_call_transition( in_states, in_symbol_token, call_transition );
            return;
        }

        std::string out_symbol_token;
        std::string move_token;
        std::string extra;

        if ( !( in >> out_symbol_token >> move_token ) ) {
            throw std::invalid_argument{
                "Invalid transition rule: " + std::string{ rule }
            };
        }

        if ( in >> extra ) {
            throw std::invalid_argument{
                "Too many tokens in transition rule: " + std::string{ rule }
            };
        }

        auto out_state{ parse_output_state_token( out_state_token ) };
        auto out_symbol{ parse_output_symbol_token( out_symbol_token ) };
        Move move{ parse_move_token( move_token ) };

        Rule_rhs rhs{ out_state, out_symbol, move };

        assign_step_transition( in_states, in_symbol_token, rhs );
    }

    void assign_step_transition(
        State_selector const &in_states,
        std::string const &in_symbol_token,
        Rule_rhs const &rhs
    ) {
        make_proxy(
            in_states,
            parse_symbol_selector( in_symbol_token ),
            false
        ) = rhs;
    }

    void assign_call_transition(
        State_selector const &in_states,
        std::string const &in_symbol_token,
        Call_transition const &call_transition
    ) {
        make_proxy(
            in_states,
            parse_symbol_selector( in_symbol_token ),
            false
        ) = call_transition;
    }

    void register_program( std::string name ) {
        if ( !std::regex_match( name, identifier_regex ) ) {
            throw std::invalid_argument{
                "'" + name + "' is not a valid program identifier."
            };
        }

        auto [it, inserted] = registry_.emplace( name, this );

        if ( !inserted ) {
            throw std::invalid_argument{
                "A Turing program named '" + name + "' is already registered."
            };
        }

        name_ = std::move( name );
    }

    static Turing_program const &find_program( std::string const &name ) {
        auto it{ registry_.find( name ) };

        if ( it == registry_.end() ) {
            throw std::invalid_argument{
                "No Turing program named '" + name + "' has been registered."
            };
        }

        return *( it->second );
    }

    inline static std::regex const identifier_regex{
        R"([A-Za-z_][A-Za-z0-9_]*)"
    };

    std::size_t state_width() const {
        std::size_t width{ 1 };

        for ( state_type s : terminal_states_ ) {
            width = std::max( width, state_to_string( s ).size() );
        }

        for ( auto const &[key, transition] : transitions_ ) {
            width = std::max( width, state_to_string( key.state ).size() );

            if ( std::holds_alternative<Step_transition>( transition ) ) {
                auto const &t{ std::get<Step_transition>( transition ) };
                width = std::max( width, state_to_string( t.new_state ).size() );
            }
        }

        return width;
    }

    void print_trace_tape(
        std::ostream &out,
        tape_type const &tape,
        typename tape_type::position_type first,
        typename tape_type::position_type last
    ) const {
        out << "~~ ";

        for ( auto pos{ first }; pos <= last; ++pos ) {
            auto code{ tape.read_code_at( pos ) };
            std::string symbol{
                code == tape_type::BLANK
                    ? "_"
                    : alphabet_.code_to_string( code )
            };

            if ( pos == tape.head() ) {
                out << '[' << symbol << ']';
            } else {
                out << ' ' << symbol << ' ';
            }

            if ( pos != last ) {
                out << ' ';
            }
        }

        out << " ~~";
    }
};

#endif
