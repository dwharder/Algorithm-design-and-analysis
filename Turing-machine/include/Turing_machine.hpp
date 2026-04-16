#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <thread>
#include <chrono>
#include <string>
#include <sstream>

class Turing_machine;
std::ostream &operator<<( std::ostream &out, Turing_machine const &tm );

class Turing_machine {
public:
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    using state_type      = size_type;
    using symbol_type     = size_type;

    static constexpr symbol_type BLANK =
        std::numeric_limits<symbol_type>::max();

    static constexpr symbol_type ANY =
        std::numeric_limits<symbol_type>::max() - 1;

    static constexpr state_type INITIAL_STATE = 0;

    static constexpr state_type HALT_STATE =
        std::numeric_limits<state_type>::max();

    enum class Move {
        left,
        right,
        stay
    };

    static constexpr Move LEFT  = Move::left;
    static constexpr Move RIGHT = Move::right;
    static constexpr Move STAY  = Move::stay;

    struct Transition {
        state_type new_state;
        symbol_type new_symbol;
        Move move;
    };

    friend std::ostream &operator<<( std::ostream &out,
                                     Turing_machine const &tm );

private:
    struct Key {
        state_type  state;
        symbol_type symbol;

        bool operator==( Key const &other ) const {
            return (state == other.state) && (symbol == other.symbol);
        }
    };

    struct Key_hash {
        size_type operator()( Key const &key ) const {
            return std::hash<state_type>{}( key.state )
                   ^ (std::hash<symbol_type>{}( key.symbol ) << 1);
        }
    };

    size_type alphabet_size_;
    state_type state_;
    difference_type head_;
    difference_type initial_head_;

    std::unordered_map<Key, Transition, Key_hash>    transitions_;
    std::unordered_map<difference_type, symbol_type> tape_;

    difference_type min_touched_{ 0 };
    difference_type max_touched_{ 0 };

    void note_access( difference_type position ) {
        if ( position < min_touched_ ) {
            min_touched_ = position;
        }

        if ( position > max_touched_ ) {
            max_touched_ = position;
        }
    }

    void validate_visible_symbol( symbol_type symbol ) const {
        if ( symbol >= alphabet_size_ ) {
            if ( (symbol == BLANK) || (symbol == ANY) ) {
                throw std::out_of_range{
                    "The BLANK symbol or ANY wildcard was used when "
                    "expecting a symbol from the visible alphabet."
                };
            } else {
                throw std::out_of_range{
                    "The symbol is outside the visible alphabet."
                };
            }
        }
    }

    void validate_transition_symbol( symbol_type symbol ) const {
        if ( symbol == ANY ) {
            throw std::out_of_range{
                "The ANY wildcard was used where a "
                "visible symbol or BLANK was required."
            };
        }

        if ( symbol == BLANK ) {
            return;
        }

        if ( symbol >= alphabet_size_ ) {
            throw std::out_of_range{
                "The symbol is outside the visible alphabet."
            };
        }
    }

    symbol_type read_raw_at( difference_type position ) {
        note_access( position );
        auto iter{ tape_.find( position ) };

        if ( iter == tape_.end() ) {
            return BLANK;
        }

        return iter->second;
    }

    symbol_type read_raw_at( difference_type position ) const {
        auto iter{ tape_.find( position ) };

        if ( iter == tape_.end() ) {
            return BLANK;
        }

        return iter->second;
    }

    void write_raw_at( difference_type position, symbol_type symbol ) {
        validate_transition_symbol( symbol );
        note_access( position );

        if ( symbol == BLANK ) {
            tape_.erase( position );
        } else {
            tape_[position] = symbol;
        }
    }

    void print_internal(
        std::ostream &out = std::cout,
        difference_type padding = 3,
        bool show_head = true
    ) const {
        difference_type first{ min_touched_ - padding };
        difference_type last{  max_touched_ + padding };

        for ( difference_type pos{ first }; pos <= last; ++pos ) {
            symbol_type symbol = read_raw_at( pos );

            if ( symbol == BLANK ) {
                out << "-";
            } else {
                out << symbol;
            }

            if ( pos != last ) {
                out << " ";
            }
        }

        if ( show_head ) {
            out << std::endl;

            std::size_t offset{ 0 };

            for ( difference_type pos{ first }; pos < head_; ++pos ) {
                symbol_type symbol = read_raw_at( pos );

                if ( symbol == BLANK ) {
                    offset += 1;
                } else {
                    offset += std::to_string( symbol ).size();
                }

                offset += 1;
            }

            for ( std::size_t k{ 0 }; k < offset; ++k ) {
                out << " ";
            }

            out << "^";
        }
    }

    state_type parse_state( std::string const &text ) const {
        if ( text == "Z" ) {
            return HALT_STATE;
        }

        std::size_t pos{ 0 };
        state_type value{};

        try {
            value = static_cast<state_type>( std::stoull( text, &pos ) );
        } catch ( ... ) {
            throw std::invalid_argument{
                "The state '" + text + "' is invalid."
            };
        }

        if ( pos != text.size() ) {
            throw std::invalid_argument{
                "The state '" + text + "' is invalid."
            };
        }

        return value;
    }

    symbol_type parse_symbol( std::string const &text ) const {
        if ( text == "*" ) {
            return ANY;
        }

        if ( text == "_" ) {
            return BLANK;
        }

        std::size_t pos{ 0 };
        symbol_type value{};

        try {
            value = static_cast<symbol_type>( std::stoull( text, &pos ) );
        } catch ( ... ) {
            throw std::invalid_argument{
                "The symbol '" + text + "' is invalid."
            };
        }

        if ( pos != text.size() ) {
            throw std::invalid_argument{
                "The symbol '" + text + "' is invalid."
            };
        }

        validate_transition_symbol( value );
        return value;
    }

    Move parse_move( std::string const &text ) const {
        if ( text == "L" ) {
            return LEFT;
        }

        if ( text == "R" ) {
            return RIGHT;
        }

        if ( text == "S" ) {
            return STAY;
        }

        throw std::invalid_argument{
            "The move '" + text + "' is invalid."
        };
    }

public:
    class Transition_proxy {
    private:
        Turing_machine &machine_;
        state_type      state_;
        symbol_type     symbol_;

    public:
        Transition_proxy(
            Turing_machine &machine,
            state_type state,
            symbol_type symbol
        ): machine_{ machine },
           state_{ state },
           symbol_{ symbol }
        {
            // Intentionally empty:
            //   no additional initialization is
            //   required beyond the initializer list.
        }

        Transition_proxy &operator=( Transition const &transition ) {
            if ( symbol_ == ANY ) {
                for (
                    symbol_type s{ 0 };
                    s < machine_.alphabet_size_;
                    ++s
                ) {
                    symbol_type out = (transition.new_symbol == ANY)
                        ? s : transition.new_symbol;

                    machine_.validate_transition_symbol( out );

                    machine_.transitions_[Key{state_, s}] = Transition{
                        transition.new_state, out, transition.move
                    };
                }

                // Handle BLANK separately
                symbol_type out = (transition.new_symbol == ANY)
                    ? BLANK : transition.new_symbol;

                machine_.validate_transition_symbol( out );

                machine_.transitions_[Key{state_, BLANK}] = Transition{
                    transition.new_state, out, transition.move
                };
            } else {
                machine_.validate_transition_symbol( symbol_ );
                machine_.validate_transition_symbol(
                    transition.new_symbol
                );

                machine_.transitions_[Key{state_, symbol_}]
                    = transition;
            }

            return *this;
        }

        operator Transition() const {
            auto iter{ machine_.transitions_.find(
                Key{state_, symbol_}
            ) };

            if ( iter == machine_.transitions_.end() ) {
                throw std::runtime_error{
                    "No transition is defined for the "
                    "specified state and symbol."
                };
            }

            return iter->second;
        }
    };

    explicit Turing_machine(
        size_type alphabet_size,
        std::initializer_list<symbol_type> initial_tape = {},
        difference_type head_position = 0
    ): alphabet_size_{ alphabet_size },
       state_{ 0 },
       head_{ head_position },
       initial_head_{ head_position },
       min_touched_{ head_position },
       max_touched_{ head_position }
    {

        if ( alphabet_size_ == 0 ) {
            throw std::invalid_argument{
                "The alphabet size must be positive."
            };
        }

        if ( alphabet_size_ >= ANY ) {
            throw std::invalid_argument{
                "The alphabet size cannot be larger than "
                    + std::to_string( ANY ) + "."

            };
        }

        difference_type pos{ 0 };

        for ( symbol_type symbol : initial_tape ) {
            validate_visible_symbol( symbol );
            write_raw_at( pos, symbol );
            ++pos;
        }
    }

    Transition_proxy operator()( state_type state, symbol_type symbol ) {
        return Transition_proxy{ *this, state, symbol };
    }

    state_type state() const {
        return state_;
    }

    bool halted() const {
        return state_ == HALT_STATE;
    }

    symbol_type read() const {
        return read_raw_at( head_ );
    }

    bool is_blank() const {
        return read_raw_at( head_ ) == BLANK;
    }

    void write( symbol_type symbol ) {
        validate_visible_symbol( symbol );
        write_raw_at( head_, symbol );
    }

    void erase() {
        write_raw_at( head_, BLANK );
    }

    void move_left() {
        --head_;
        note_access( head_ );
    }

    void move_right() {
        ++head_;
        note_access( head_ );
    }

    void set_state( state_type new_state ) {
        state_ = new_state;
    }

    void step() {
        if ( halted() ) {
            throw std::logic_error{
                 "Cannot step a machine that has already halted."
            };
        }

        symbol_type current_symbol = read_raw_at( head_ );
        auto iter{ transitions_.find( Key{state_, current_symbol} ) };

        if ( iter == transitions_.end() ) {
            throw std::runtime_error{
                "No transition is defined for the "
                "current state and symbol."
            };
        }

        Transition const &transition = iter->second;

        write_raw_at( head_, transition.new_symbol );
        state_ = transition.new_state;

        switch ( transition.move ) {
        case Move::left:
            move_left();
            break;
        case Move::right:
            move_right();
            break;
        case Move::stay:
            break;
        }
    }

    void run( size_type max_steps = 100000 ) {
        for (
            size_type count{ 0 };
            (count < max_steps) && !halted();
            ++count
        ) {
            step();
        }

        if ( !halted() ) {
            throw std::runtime_error{
                "The machine did not halt within "
                "the maximum number of steps."
            };
        }
    }

    void run_with_delay(
        int delay_ms = 200,
        std::ostream &out = std::cout,
        difference_type padding = 3,
        bool show_head = true,
        size_type max_steps = 100000
    ) {
        print_internal( out, padding, show_head );
        out << std::endl;

        for (
            size_type count{ 0 };
            (count < max_steps) && !halted();
            ++count
        ) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds( delay_ms )
            );

            step();

            print_internal( out, padding, show_head );
            out << std::endl;
        }

        if ( !halted() ) {
            throw std::runtime_error(
                "The machine did not halt within "
                "the maximum number of steps."
            );
        }
    }

    void run_from(
        state_type start_state,
        size_type max_steps = 100000
    ) {
        state_ = start_state;
        run( max_steps );
    }

    void reset() {
        state_ = INITIAL_STATE;
        head_  = initial_head_;
    }

    void print(
        std::ostream &out = std::cout,
        difference_type padding = 3,
        bool show_head = true
    ) const {
        print_internal( out, padding, show_head );
        out << std::endl;
    }

    void single_step_and_print(
        std::ostream &out = std::cout,
        difference_type padding = 3,
        bool show_head = true
    ) {
        step();
        print_internal( out, padding, show_head );
        out << std::endl;
    }

    void add_rule( std::string const &rule ) {
        std::istringstream in{ rule };

        std::string current_state_str;
        std::string current_symbol_str;
        std::string arrow;
        std::string new_state_str;
        std::string new_symbol_str;
        std::string move_str;

        if ( !(in >> current_state_str
                  >> current_symbol_str
                  >> arrow
                  >> new_state_str
                  >> new_symbol_str
                  >> move_str) ) {
            throw std::invalid_argument{
                "The rule could not be parsed."
            };
        }

        if ( arrow != "->" ) {
            throw std::invalid_argument{
                "The rule must contain '->'."
            };
        }

        std::string extra;
        if ( in >> extra ) {
            throw std::invalid_argument{
                "The rule contains extra text."
            };
        }

        state_type   current_state{ parse_state( current_state_str ) };
        symbol_type current_symbol{ parse_symbol( current_symbol_str ) };
        state_type       new_state{ parse_state( new_state_str ) };
        symbol_type     new_symbol{ parse_symbol( new_symbol_str ) };
        Move                  move{ parse_move( move_str ) };

        (*this)( current_state, current_symbol ) = {
            new_state, new_symbol, move
        };
    }

    void add_rules( std::initializer_list<std::string> rules ) {
        for ( std::string const &rule : rules ) {
            add_rule( rule );
        }
    }
};

std::ostream &operator<<( std::ostream &out, Turing_machine const &tm ) {
    tm.print_internal( out, 3, true );
    return out;
}
