#include <cstddef>      // std::size_t
#include <limits>       // std::numeric_limits
#include <type_traits>  // std::is_integral_v, std::is_unsigned_v
#include <vector>       // std::vector
#include <cassert>

template <typename T>
std::vector<T> primes_up_to( T n );

/**
 * @file   sieve.hpp
 * @brief  Implements the Sieve of Eratosthenes algorithm.
 *
 * This header provides a templated function to compute all prime numbers up to a given limit using the Sieve of Eratosthenes.
 * It requires the template type to be an unsigned integral type.
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
  std::vector<bool> is_prime( N + 1U, true );
  is_prime[0] = false;
  is_prime[1] = false;

  std::size_t n_primes{ 0 };

  // Outer loop condition written as k <= N/k to avoid k*k overflow
  std::size_t k{ 2U };

  for ( ; k <= N / k; ++k ) {
    if ( is_prime[k] ) {
      ++n_primes;

      // Start at k^2; advance by repeated addition (m += k).
      for ( std::size_t m = k*k; m <= N; m += k ) {
        is_prime[m] = false;
      }
    }
  }

  // Count the remaining prime numbers
  for ( ; k <= N; ++k ) {
    if ( is_prime[k] ) {
      ++n_primes;
    }
  }

  std::vector<T> primes( n_primes );
  std::size_t m{ 0 };

  for ( std::size_t k{ 2U }; m < n_primes; ++k ) {
    assert( k <= N );

    if ( is_prime[k] ) {
      primes[m++] = static_cast<T>( k );
    }
  }

  return primes;
}
