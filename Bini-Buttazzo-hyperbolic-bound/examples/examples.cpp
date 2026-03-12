#include <iostream>
#include "hyperbolic_bound.hpp"


int main() {
    std::cout << std::hex;

    for ( unsigned int k{ 1 }; k <= 10; ++k ) {
        std::cout << Hyperbolic_bound::encode_utilization_factor( k, 11 ) << std::endl;
    }

    // A quick test: verified with Maple using rational arithmetic
    // and appropriate rounding, the product
    //
    //     \prod_{k = 1}^{77162} (1 + k2^{-32}) = 1.11111111111111111101101010100100
    //                                          = 1.ffffdaa4

    Hyperbolic_bound util_check{ 1 };

    std::uint32_t k{2};

    for ( ; util_check.add_task( k ); ++k );

    std::cout << std::dec << k << std::endl;
    std::cout << std::hex << util_check.raw_encoding() << std::endl;

    std::cout << std::endl;

    Hyperbolic_bound remove_check{ 0xff358732u };
    std::cout << std::hex << remove_check.raw_encoding() << std::endl;
    remove_check.remove_task( 0x159ef900u );
    std::cout << std::hex << remove_check.raw_encoding() << std::endl;
    auto result{ remove_check.add_task( 0x159ef900u ) };
    std::cout << std::hex << remove_check.raw_encoding() << std::endl;

    return 0;
}
