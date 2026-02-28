#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <array>

// In a partitioning algorithm (e.g., quicksort or quickselect),
// the optimal pivot is the true median of the range, since this
// produces two equal-sized partitions. However, computing the
// median exactly requires an O(n) operation, which defeats the
// purpose of fast divide-and-conquer schemes. In practice, one
// instead approximates the median by sampling a small subset
// of the range.
//
// This program compares four pivot-selection strategies:
//
//   1. A single randomly chosen element.
//   2. The larger of two randomly chosen elements.
//   3. The median of three randomly chosen elements.
//   4. Tukey's ninther: the median of three medians, each formed
//      from three randomly chosen elements (nine samples total).
//
// The performance of divide-and-conquer algorithms depends
// primarily on the size of the larger partition. Although the
// expected behavior of these strategies can be derived
// analytically, this program performs a Monte Carlo simulation
// to confirm that empirical results align with the theoretical
// predictions.
//
// For each trial, a vector of N random values is generated and
// partitioned using each strategy. To ensure a fair comparison,
// all strategies draw from the same sampled indices whenever
// their definitions overlap. This prevents any strategy from
// benefiting merely from a fortunate random selection within
// a given trial.

  ///////////////////////////
 // Function declarations //
///////////////////////////
int main();

// Return the size of the larger side after partitioning
// around pivot_value, using predicate x < pivot_value
//  - so "left" is < pivot, "right" is >= pivot
template <class T>
std::size_t larger_side_after_partition(
    std::vector<T> &v, T const &pivot_value
);

// Median of three values using only operator< (no ==, no <=).
template <class T>
T const &median_of_three( T const &a, T const &b, T const &c );

  ////////////////////////////
 // Function dedefinitions //
////////////////////////////

int main() {
    // Parameters:
    constexpr std::size_t N{ 10000 };       // vector size
    constexpr std::size_t TRIALS{ 2000 };   // number of trials

    // Random engine:
    std::mt19937_64 rng{ std::random_device{}() };

    // Unsigned integer distribution (full range)
    std::uniform_int_distribution<unsigned int> value_dist(
        0u, std::numeric_limits<unsigned int>::max()
    );

    // Index distribution (0..N-1)
    std::uniform_int_distribution<std::size_t> index_dist( 0, N - 1 );

    double sum_larger_one_pivot{ 0.0 };
    double sum_larger_two_pivot{ 0.0 };
    double sum_median_of_three{ 0.0 };
    double sum_median_of_medians{ 0.0 };

    for ( std::size_t trial{ 0 }; trial < TRIALS; ++trial ) {
        // Generate fresh random data each trial
        //  - Avoid unnecessarily initializing 'a'
        std::vector<unsigned int> a;
        a.reserve( N );

        for ( std::size_t i{ 0 }; i < N; ++i ) {
            a.push_back( value_dist( rng ) );
        }

        // Make copies to be partitioned
        std::vector<unsigned int> v1{ a };
        std::vector<unsigned int> v2{ a };
        std::vector<unsigned int> v3{ a };
        std::vector<unsigned int> v4{ a };

        // So that each strategy uses the same indices, so that one
        // strategy does not get an advantage over another by random chance.
        //  - Weakness, we don't check for duplicates
        std::array<std::size_t, 9> indices;
      
        for ( auto &idx : indices ) {
            idx = index_dist( rng );
        }

        // Strategy 1: pivot is one randomly chosen element
        {
            unsigned int pivot{ v1[indices[0]] };

            std::size_t larger{ larger_side_after_partition( v1, pivot ) };
            sum_larger_one_pivot += static_cast<double>( larger );
        }

        // Strategy 2: pivot is the larger of two randomly chosen elements
        {
            unsigned int p1{ v2[indices[0]] };
            unsigned int p2{ v2[indices[1]] };
            unsigned int pivot{ (p1 < p2) ? p2 : p1 };  // max(p1,p2)

            std::size_t larger{ larger_side_after_partition( v2, pivot ) };
            sum_larger_two_pivot += static_cast<double>( larger );
        }

        // Strategy 3: pivot is the median of three randomly chosen elements
        {
            unsigned int const &p1{ v3[indices[0]] };
            unsigned int const &p2{ v3[indices[1]] };
            unsigned int const &p3{ v3[indices[2]] };

            unsigned int pivot{ median_of_three( p1, p2, p3 ) };
            std::size_t larger{ larger_side_after_partition( v3, pivot ) };
            sum_median_of_three += static_cast<double>( larger );
        }

        // Strategy 4: pivot is the median of three medians of three
        // randomly chosen elements (Tukey's ninther)
        {
            // Pick 9 samples (with replacement, as in the model).
            unsigned int const &a1{ v4[indices[0]] };
            unsigned int const &a2{ v4[indices[1]] };
            unsigned int const &a3{ v4[indices[2]] };

            unsigned int const &b1{ v4[indices[3]] };
            unsigned int const &b2{ v4[indices[4]] };
            unsigned int const &b3{ v4[indices[5]] };

            unsigned int const &c1{ v4[indices[6]] };
            unsigned int const &c2{ v4[indices[7]] };
            unsigned int const &c3{ v4[indices[8]] };

            unsigned int m1{ median_of_three( a1, a2, a3 ) };
            unsigned int m2{ median_of_three( b1, b2, b3 ) };
            unsigned int m3{ median_of_three( c1, c2, c3 ) };

            unsigned int pivot{ median_of_three( m1, m2, m3 ) };

            std::size_t larger{ larger_side_after_partition( v4, pivot ) };
            sum_median_of_medians += static_cast<double>( larger );
        }
    }

    double avg_one{  sum_larger_one_pivot/static_cast<double>( TRIALS ) };
    double avg_two{  sum_larger_two_pivot/static_cast<double>( TRIALS ) };
    double avg_med{   sum_median_of_three/static_cast<double>( TRIALS ) };
    double avg_m3m{ sum_median_of_medians/static_cast<double>( TRIALS ) };

    std::cout << "Size: " << N << std::endl;
    std::cout << "Number of trials: " << TRIALS << std::endl;
    std::cout << std::fixed << std::setprecision( 3 );

    std::cout << "Average larger side (one random pivot):        "
              << avg_one
              << "  (expected ratio: 0.75,   actual ratio: "
              << (avg_one / static_cast<double>( N ))
              << ")" << std::endl;

    std::cout << "Average larger side (larger of two pivots):    "
              << avg_two
              << "  (expected ratio: 0.75,   actual ratio: "
              << (avg_two / static_cast<double>( N ))
              << ")" << std::endl;

    std::cout << "Average larger side (median of three):         "
              << avg_med
              << "  (expected ratio: 0.6875, actual ratio: "
              << (avg_med / static_cast<double>( N ))
              << ")" << std::endl;

    std::cout << "Average larger side (median of three medians): "
              << avg_m3m
              << "  (expected ratio: 0.6336, actual ratio: "
              << (avg_m3m / static_cast<double>( N ))
              << ")" << std::endl;

    return 0;
}

// Return the size of the larger side after partitioning
// around pivot_value, using predicate x < pivot_value
//  - So "left" is < pivot, "right" is >= pivot
template <class T>
std::size_t larger_side_after_partition(
    std::vector<T> &v, T const &pivot_value
) {
    auto it = std::partition( v.begin(), v.end(), [&]( T const &x ) {
        return x < pivot_value;
    } );

    std::size_t left{
        static_cast<std::size_t>( std::distance( v.begin(), it ) )
    };

    std::size_t right{ v.size() - left };
    return ( left > right ) ? left : right;
}

// Median of three values using only operator< (no ==, no <=).
template <class T>
T const &median_of_three( T const &a, T const &b, T const &c ) {
    // Returns the middle element of {a,b,c}.
    if ( a < b ) {
        if ( b < c ) {
            return b;         // a < b < c
        } else if ( a < c ) {
            return c;         // a < c <= b
        } else {
            return a;         // c <= a < b
        }
    } else { // b <= a
        if ( a < c ) {
            return a;         // b <= a < c
        } else if ( b < c ) {
            return c;         // b < c <= a
        } else {
            return b;         // c <= b <= a
        }
    }
}
