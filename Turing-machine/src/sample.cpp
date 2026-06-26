#include <iostream>
#include "Turing.hpp"

using namespace TM;

int main() {
    Alphabet Binary{ '0', '1'  };

    Turing_program Move_to_left_most_symbol{ "move_to_left_most_symbol", Binary, {
        "0 _ -> 0 _ L",   // Keep moving to the left until you reach a symbol
        "0 # -> 1 = L",   // Having reached a symbol, move to the left and -> State 1
        "1 # -> 1 = L",   // So long as you are over a symbol, move to the left
        "1 _ -> Z _ R"    // Having reached a blank, move right and halt
    } };

    Turing_program Move_to_right_most_symbol{ "move_to_right_most_symbol", Binary, {
        "0 _ -> 0 _ R",   // Keep moving to the right until you reach a symbol
        "0 # -> 1 = R",   // Having reached a symbol, move to the right and -> State 1
        "1 # -> 1 = R",   // So long as you are over a symbol, move to the right
        "1 _ -> Z _ L"    // Having reached a blank, move left and halt
    } };

    Turing_program Is_non_negative_binary{ "Is_non_negative_binary", Binary, {
        "0   _   -> F _ S",
        "0 {0,1} -> call move_to_right_most_symbol Z->1",
        "1 0 -> 2 0 L",
        "1 1 -> 4 1 L",
        "2 _ -> T _ R",
        "2 0 -> 3 0 L",
        "2 1 -> 4 1 L",
        "3 _ -> 5 _ R",
        "3 0 -> 3 0 L",
        "3 1 -> 4 1 L",
        "4 _ -> 6 _ R",
        "4 0 -> 3 0 L",
        "4 1 -> 4 1 L",
        "5 {0,1} -> call move_to_right_most_symbol Z->F",
        "6 {0,1} -> call move_to_right_most_symbol Z->T"
    }, {TRUE, FALSE} };

    Turing_program Is_byte{ "Is_byte", Binary, {
        "0   _   -> F _ S",
        "0 {0,1} -> call move_to_left_most_symbol Z->1",
        "1 {0,1} -> 2 = R",
        "1   _   -> F = L",
        "2 {0,1} -> 3 = R",
        "2   _   -> F = L",
        "   3 {0,1} -> 4 = R    ",
        "3   _   -> F = L",
        "4 {0,1} -> 5 = R",
        "4   _   -> F = L",
        "5 {0,1} -> 6 = R",
        "5   _   -> F = L",
        "6 {0,1} -> 7 = R",
        "6   _   -> F = L",
        "7 {0,1} -> 8 = R",
        "7   _   -> F = L",
        "8 {0,1} -> 9 = R",
        "8   _   -> F = L",
        "9   _   -> T = L",
        "9 {0,1} -> call move_to_right_most_symbol Z->F"
    }, {TRUE, FALSE} };

    Turing_program Complement{ "Complement", Binary, {
        "0 {0,1} -> call move_to_left_most_symbol Z->1",
        "1 0 -> 1 1 R",
        "1 1 -> 1 0 R",
        "1 _ -> Z _ L"
    } };

    Turing_program Twos_complement{ "Twos_complement", Binary, {
        "0 {0,1} -> call move_to_right_most_symbol Z->1",
        "1 0 -> 1 0 L",
        "1 1 -> 2 1 L",
        "1 _ -> call move_to_right_most_symbol",
        "2 0 -> 2 1 L",
        "2 1 -> 2 0 L",
        "2 _ -> call move_to_right_most_symbol",
    } };

    Tape t1{ Binary, "01011100", 3 }, t2{ Binary, "11011100", 4 }, t3{ Binary, "11011101", 5 };
    std::cout << t1 << " " << t2 << " " << t3 << std::endl;

    std::cout << "All of these are bytes:" << std::endl;
    std::cout << Is_byte( t1 ) << std::endl;
    std::cout << Is_byte( t2 ) << std::endl;
    std::cout << Is_byte( t3 ) << std::endl << std::endl;

    std::cout << "The first does not have a leading '1' and is not zero:" << std::endl;
    std::cout << Is_non_negative_binary( t1 ) << std::endl;
    std::cout << Is_non_negative_binary( t2 ) << std::endl;
    std::cout << Is_non_negative_binary( t3 ) << std::endl << std::endl;

    std::cout << "The complement and 2's complement of these bytes:" << std::endl;
    std::cout << t1 << std::endl;
    Complement( t1 );
    std::cout << t1 << std::endl;
    Complement( t1 );
    Twos_complement( t1 );
    std::cout << t1 << std::endl << std::endl;

    std::cout << t2 << std::endl;
    Complement( t2 );
    std::cout << t2 << std::endl;
    Complement( t2 );
    Twos_complement( t2 );
    std::cout << t2 << std::endl << std::endl;

    std::cout << t3 << std::endl;
    Complement( t3 );
    std::cout << t3 << std::endl;
    Complement( t3 );
    Twos_complement( t3 );
    std::cout << t3 << std::endl << std::endl;

    Tape t4{ Binary,  "0" }, t5{ Binary,  "1" }, t6{ Binary, "10" },
         t7{ Binary, "01" }, t8{ Binary, "11" };
    std::cout << t4 << " " << t5 << " " << t6 << " " << t7 << " " << t8 << std::endl;

    std::cout << "None of these are bytes:" << std::endl;
    std::cout << Is_byte( t4 ) << std::endl;
    std::cout << Is_byte( t5 ) << std::endl;
    std::cout << Is_byte( t6 ) << std::endl;
    std::cout << Is_byte( t7 ) << std::endl;
    std::cout << Is_byte( t8 ) << std::endl;

    std::cout << "The fourth does not have a leading '1' and is not zero:" << std::endl;
    std::cout << Is_non_negative_binary( t4 ) << std::endl;
    std::cout << Is_non_negative_binary( t5 ) << std::endl;
    std::cout << Is_non_negative_binary( t6 ) << std::endl;
    std::cout << Is_non_negative_binary( t7 ) << std::endl;
    std::cout << Is_non_negative_binary( t8 ) << std::endl;

    return 0;
}
