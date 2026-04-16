# Turing Machine Library (C++)

This library provides a simple, structured interface for constructing and executing Turing machines in C++.

It supports three approaches:

1. **Transition rules using proxy assignments** (recommended)  
2. **A domain-specific language (DSL)** for compact rule definitions  
3. **Manual control** of the machine (for experimentation and debugging)  

---

## 1. Transition Rules (Proxy Assignment)

This approach is closest to the formal definition of a Turing machine.

### Example

    #include "Turing_machine.hpp"

    int main() {
        using TM = Turing_machine;

        TM tm{ 2, {0,0,1,0,1,0,1,0}, 7 };

        // Complement each bit while moving left
        for ( std::size_t s{ 0 }; s < 8; ++s ) {
            tm( s, 0 ) = { s + 1, 1, TM::LEFT };
            tm( s, 1 ) = { s + 1, 0, TM::LEFT };
        }

        // Move right 8 times
        //  - Any symbol or blank should be left unchanged
        for ( std::size_t s{ 8 }; s < 16; ++s ) {
            tm( s, TM::ANY ) = { s + 1, TM::ANY, TM::RIGHT };
        }

        // Halt
        //  - Go to the halting state leaving the symbol or blank unchanged
        tm( 16, TM::ANY ) = { TM::HALT_STATE, TM::ANY, TM::STAY };

        tm.run();
    }

### Transition Definition

    tm( state, symbol ) = { new_state, new_symbol, move };

This corresponds to:

    δ(q, a) = (q', b, d)

### Constants

- `TM::BLANK` — the blank symbol  
- `TM::ANY` — wildcard (matches any symbol)  
- `TM::HALT_STATE` — halting state  
- `TM::LEFT`, `TM::RIGHT`, `TM::STAY` — head movement  

---

## 2. DSL (Domain-Specific Language)

A compact string-based format for defining transitions.

### Example

    Turing_machine tm{ 2, {0,0,1,0,1,0,1,0}, 7 };

    tm.add_rules({
        " 0 0 ->  1 1 L",
        " 0 1 ->  1 0 L",
        " 1 0 ->  2 1 L",
        " 1 1 ->  2 0 L",
        " 2 0 ->  3 1 L",
        " 2 1 ->  3 0 L",
        " 3 0 ->  4 1 L",
        " 3 1 ->  4 0 L",
        " 4 0 ->  5 1 L",
        " 4 1 ->  5 0 L",
        " 5 0 ->  6 1 L",
        " 5 1 ->  6 0 L",
        " 6 0 ->  7 1 L",
        " 6 1 ->  7 0 L",
        " 7 0 ->  8 1 L",
        " 7 1 ->  8 0 L",

        " 8 * ->  9 * R",
        " 9 * -> 10 * R",
        "10 * -> 11 * R",
        "11 * -> 12 * R",
        "12 * -> 13 * R",
        "13 * -> 14 * R",
        "14 * -> 15 * R",
        "15 * -> 16 * R",

        "16 * ->  Z * S"
    });

    tm.run();

### DSL Syntax

    state symbol -> new_state new_symbol move

### Special Symbols

- `*` — wildcard (`TM::ANY`)  
- `_` — blank (`TM::BLANK`)  
- `Z` — halting state (`TM::HALT_STATE`)  
- `L`, `R`, `S` — movement  

### Notes

- `*` matches all symbols (including blank)  
- `*` on the right-hand side preserves the symbol  
- Later rules overwrite earlier ones  

---

## 3. Manual Control

You may directly manipulate the machine:

    Turing_machine tm{ 2, {0,1,1,0}, 3 };

    auto symbol = tm.read();

    if ( symbol == 0 ) {
        tm.write( 1 );
    } else {
        tm.write( 0 );
    }

    tm.move_left();
    tm.erase();

### Available Operations

- `read()` — returns current symbol (may be `BLANK`)  
- `write(symbol)` — write visible symbol  
- `erase()` — write blank  
- `move_left()`, `move_right()`  
- `step()` — execute one transition  
- `run()` — execute until halt  
- `run_with_delay(ms)` — step-by-step execution  
- `reset()` — reset state and head (not tape)  

---

## Output

    - - - 0 0 1 0 1 0 1 0 - - -
                        ^

- `-` represents blank  
- `^` marks the head position  

---

## Design Notes

- The tape is conceptually infinite and stored sparsely  
- Transitions are stored using hash maps  
- Proxy assignments provide a clean interface  
- `ANY` expands to all symbols, including blank  

---

## Acknowledgement

The overall design of this package, including the use of proxy transitions, transition tables based on unordered maps, and support for wildcards such as `TM::ANY`, was developed by the author.

The implementation was developed collaboratively with ChatGPT, which contributed significantly to the code, suggested refinements, and helped develop features such as the DSL interface.
