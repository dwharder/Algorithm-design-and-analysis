#include <cstddef>      // std::size_t
#include <limits>       // std::numeric_limits
#include <type_traits>  // std::is_integral_v, std::is_unsigned_v
#include <vector>       // std::vector
#include <cassert>

// Forward declaration to allow inclusion as a header-only implementation.
template <typename T>
std::vector<T> primes_up_to( T n );

/**
 * @file   sieve.hpp
 * @brief  Implements the Sieve of Eratosthenes algorithm.
 *
 * This header provides a templated function to compute all prime
 * numbers up to a given limit using the Sieve of Eratosthenes.
 *     Time complexity:  O(n ln(ln(n)) )
 *     Space complexity: O(n)
 * \link https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes the Wikipedia article \endlink.
 *
 * This requires the template type to be an unsigned integral type.
 * This header is intended to be header-only and contains the full implementation.
 *
 * @tparam T An unsigned integral type representing the upper limit of primes to find.
 * @param  n The upper limit up to which primes are computed.
 * @return A std::vector<T> containing all prime numbers up to n.
 */

template <typename T>
std::vector<T> primes_up_to( T n ) {
  static_assert( std::is_integral_v<T>,
                 "primes_up_to<T>: T must be an integral type." );

  static_assert( std::is_unsigned_v<T>,
                 "primes_up_to<T>: T must be an unsigned integral type." );

  // No primes below 2
  if ( n < T{ 2 } ) {
    return {};
  }

  const std::size_t N{ static_cast<std::size_t>( n ) };

  // is_prime[k] == true means "k is prime" (until proven composite)
  // We use std::vector<bool> for compact storage
  //  - bit-level representation is acceptable here
  //    since we only require indexed reads and writes.
  std::vector<bool> is_prime( N + 1, true );
  is_prime[0] = false;
  is_prime[1] = false;

  // Total number of primes <= n
  //  - Computed first to allow exact allocation.
  std::size_t n_primes{ 0 };

  std::size_t k{ 2 };

  // Outer loop condition written as k <= N/k to avoid k*k overflow
  //  - Equivalent to k*k <= N for positive integers
  for ( ; k <= N / k; ++k ) {
    if ( is_prime[k] ) {
      ++n_primes;

      // Start at k*k: advance by repeated addition (m += k).
      for ( std::size_t m = k*k; m <= N; m += k ) {
        is_prime[m] = false;
      }
    }
  }

  // Count primes greater than sqrt(n) that were not visited in the sieve loop
  for ( ; k <= N; ++k ) {
    if ( is_prime[k] ) {
      ++n_primes;
    }
  }

  #ifndef NDEBUG
    std::size_t check{ 0 };
  
    for ( std::size_t i{ 2 }; i <= N; ++i ) {
      if ( is_prime[i] ) ++check;
    }
  
    assert( check == n_primes );
  #endif

  // Allocate result vector with exact size to avoid reallocation
  std::vector<T> primes( n_primes );

  for ( std::size_t m{ 0 }, k{ 2U }; m < n_primes; ++k ) {
    // Defensive check: logic guarantees termination before k exceeds N
    assert( k <= N );

    if ( is_prime[k] ) {
      primes[m++] = static_cast<T>( k );
    }
  }

  return primes;
}
