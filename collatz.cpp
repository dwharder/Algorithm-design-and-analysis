#include <cstddef>
#include <iostream>
#include <vector>
#include <cassert>
#include <limits>
#include <iomanip>

/**
 * @file collatz.cpp
 *
 * @brief Compute record total stopping times for the Collatz map.
 *
 * @author Douglas Wilhelm Harder
 *
 * @date 2026
 *
 * @details
 * For a positive integer \(n\), the Collatz map is
 * \f[   C(n) =
 *   \begin{cases}
 *     n/2,   & \text{if } n \text{ is even}, \\
 *     3n+1,  & \text{if } n \text{ is odd}.
 *   \end{cases}   \f]
 *
 * The total stopping time \(T(n)\) is the number of iterations required
 * for repeated application of \(C\) to reach \(1\).
 *
 * This program computes \(T(n)\) for all \(1 \le n < N\), printing each
 * \(n\) for which \(T(n)\) exceeds all previously observed stopping times.
 *
 * Because intermediate values in a Collatz trajectory may exceed the range
 * of a single `std::size_t`, this program represents intermediate values
 * using two words, thereby emulating an unsigned integer of width
 * \(2 \cdot 8 \cdot \texttt{sizeof(std::size_t)}\) bits.
 *
 * @section algorithm_overview Algorithm overview
 *
 * The implementation uses dynamic programming:
 *
 * - previously computed stopping times are cached in an array `v`
 * - when the trajectory of a new value reaches a cached value, the
 *   remaining stopping time is known immediately
 * - the program then retraces the trajectory from the starting value and
 *   backfills the cache for intermediate values below `N`
 *
 * @section optimization_notes Optimizations
 *
 * The implementation uses the following optimizations:
 *
 * - if \(k\) is even, then \(T(k)=T(k/2)+1\), so no trajectory walk is needed
 * - for odd intermediate values, the two successive steps
 *   \f[   n \mapsto 3n+1 \mapsto (3n+1)/2   \f]
 *   are combined into one function call
 * - fast one-word paths are used when the intermediate value fits in a
 *   single `std::size_t`
 */

/**
 * @brief Two-word unsigned integer.
 *
 * @details
 * This structure represents the value
 * \f[   n = \texttt{high} \cdot 2^w + \texttt{low},   \f]
 * where
 * \f[   w = 8 \cdot \texttt{sizeof(std::size_t)}.   \f]
 *
 * Thus `high` stores the most significant word and `low` stores the least
 * significant word.
 *
 * @note
 * This is not a general arbitrary-precision type. It is only the minimal
 * two-word representation needed by this program. It is only when N is on
 * the order of 20 billion intermediate values are as large as
 * 0x14da9215f6616e8f0 > 0xffffffffffffffff = 2^64 - 1.
 */
struct double_size_t {
    std::size_t high;  /**< Most significant word. */
    std::size_t low;   /**< Least significant word. */
};

/**
 * @brief Mask whose only set bit is the most significant bit of `std::size_t`.
 *
 * @details
 * This value is used when shifting a two-word integer right by one bit:
 * if the least significant bit of the high word is 1, that bit must become
 * the most significant bit of the low word.
 */
static std::size_t const MSB{
    std::size_t{ 1 } << (8*sizeof(std::size_t) - 1)
};

int main();

/**
 * @brief Divide a two-word integer by two.
 *
 * @param[in,out] n
 * The value to be updated in place.
 *
 * @details
 * This function performs a logical right shift by one bit on the two-word
 * integer represented by `n`.
 *
 * If the least significant bit of `n.high` is 1 before the shift, that bit
 * becomes the most significant bit of `n.low` after the shift.
 *
 * @post
 * On return, `n` has been replaced by \f$\lfloor n/2 \rfloor\f$.
 */
void divide_n_by_two( double_size_t &n );

/**
 * @brief Replace a two-word integer by \(3n+1\).
 *
 * @param[in,out] n
 * The value to be updated in place.
 *
 * @details
 * If `n.high == 0` and the computation cannot overflow a single word,
 * a fast one-word path is used:
 * \f[   n \leftarrow 3n + 1.   \f]
 *
 * Otherwise, the result is computed using explicit two-word arithmetic.
 *
 * @post
 * On return, `n` has been replaced by \(3n+1\).
 */
void three_n_plus_one( double_size_t &n );

/**
 * @brief Replace a two-word integer by \((3n+1)/2\).
 *
 * @param[in,out] n
 * The value to be updated in place.
 *
 * @details
 * This function is intended for use only when `n` is odd.
 * For odd \(n\), the Collatz step \(3n+1\) is always even, so the two-step
 * transition
 * \f[   n \mapsto 3n+1 \mapsto \frac{3n+1}{2}   \f]
 * can be combined into one operation.
 *
 * If `n.high == 0` and the result still fits in one word, a fast path is used:
 * \f[   \frac{3n+1}{2} = n + \frac{n+1}{2}.   \f]
 *
 * Otherwise, the result is computed using explicit two-word arithmetic.
 *
 * @pre
 * `n` should be odd when this function is called.
 *
 * @post
 * On return, `n` has been replaced by \((3n+1)/2\).
 */
void three_n_plus_one_and_divide_by_2( double_size_t &n );

/**
 * @brief Program entry point.
 *
 * @return `0` on success.
 *
 * @details
 * Let \(T(n)\) denote the total stopping time of \(n\).
 * The program computes \(T(n)\) for all \(1 \le n < N\), printing each
 * value of \(n\) that sets a new record.
 *
 * The cache array `v` stores stopping times for values below `N`.
 * Its meaning is:
 *
 * - `v[x] == 0` means “unknown”, except for the special value `v[1] == 0`
 * - for all other cached entries, `v[x]` equals \(T(x)\)
 *
 * The algorithm for an odd starting value \(k\) is:
 *
 * 1. follow the Collatz trajectory until a cached value is reached
 * 2. count how many steps were taken
 * 3. restart from \(k\)
 * 4. retrace the same trajectory and fill in cache values for intermediate
 *    states below `N`
 *
 * This avoids storing the full trajectory explicitly.
 *
 * @warning
 * The declaration
 * `std::vector<unsigned short> v(N, 0);`
 * requires a very large amount of memory when `N` is large. For example,
 * \(N = 35{,}700{,}000{,}000\) would be far beyond typical RAM capacity.
 */

/**
 * Sample output (running time of approximately 15 min with N = 35.7 billion):
@verbatim
      N = 35700000000
      <snip />
      11200681        688
      14934241        691
      15733191        704
      31466382        705
      36791535        744
      63728127        949
      127456254       950
      169941673       953
      226588897       956
      268549803       964
      537099606       965
      670617279       986
      1341234558      987
      1412987847      1000
      1674652263      1008
      2610744987      1050
      4578853915      1087
      4890328815      1131
      9780657630      1132
      12212032815     1153
      12235060455     1184
      13371194527     1210
      17828259369     1213
      31694683323     1219
      Maximum n:    0x164d6de9cfb4768b0
       where 2^64 = 0x10000000000000000
@endverbatim
 */

int main() {
    // std::size_t const N{ 2000000000 };
    std::size_t const N{ 35700000000 };

    /**
     * @brief Cache of stopping times.
     *
     * @details
     * For values \(x < N\), `v[x]` stores \(T(x)\) once known.
     *
     * The type `unsigned short` is sufficient here because all stopping times
     * of interest are far below 65535.
     */
    std::vector<unsigned short> v( N, 0 );
    std::cout << "N = " << N << std::endl;
  
    v[0] = 0;
    v[1] = 0;
    v[2] = 1;
    v[3] = 7;
    v[4] = 2;

    /**
     * @brief Largest intermediate value encountered.
     *
     * @details
     * This tracks the largest value reached by any trajectory examined during
     * the computation, represented in two-word form.
     */
    double_size_t max_n{ 0, 0 };

    std::cout << "1\t0" << std::endl;
    std::cout << "2\t1" << std::endl;
    std::cout << "3\t7" << std::endl;

    /**
     * @brief Current record stopping time.
     *
     * @details
     * Whenever a new value `v[k]` exceeds `curr_max`, the pair
     * `(k, v[k])` is printed.
     */
    std::size_t curr_max{ 7 };

    for ( std::size_t k{ 5 }; k < N; ++k ) {
        /// Skip values already backfilled during an earlier trajectory.
        if ( v[k] != 0 ) {
            continue;
        }

        /**
         * @details
         * If \(k\) is even, then
         * \f[   T(k) = T(k/2) + 1,   \f]
         * and the value is obtained immediately from the cache.
         */
        if ( (k & std::size_t{ 1 }) == 0 ) {
            v[k] = v[k >> 1] + 1;
        } else {
            double_size_t n{ 0, k };

            /**
             * @brief Number of Collatz steps taken before reaching a cached value.
             */
            unsigned short count{ 0 };

            /**
             * Walk forward until the current value is cacheable and already known.
             */
            while ( (n.high != 0) || (n.low >= N) || (v[n.low] == 0) ) {
                if ( (n.low & 1u) == 0 ) {
                    divide_n_by_two( n );
                    ++count;
                } else {
                    /**
                     * For odd \(n\), compute
                     * \f[   n \leftarrow \frac{3n+1}{2}   \f]
                     * and count this as two Collatz steps.
                     */
                    three_n_plus_one_and_divide_by_2( n );
                    count += 2;
                }

                /// Update the largest intermediate value seen so far.
                if ( (n.high > max_n.high) ||
                     ( (n.high == max_n.high) && (n.low > max_n.low) ) ) {
                    max_n.high = n.high;
                    max_n.low  = n.low;
                }
            }

            /**
             * @brief Cached tail length.
             *
             * @details
             * Once the trajectory reaches a cached value `n.low`,
             * the remaining number of steps to 1 is `v[n.low]`.
             */
            unsigned short base{ v[n.low] };

            /// Restart from the original starting value.
            n.high = 0;
            n.low  = k;

            /**
             * @brief Stopping time of the current state during backfill.
             *
             * @details
             * Initially this equals `count + base`, which is exactly \(T(k)\).
             * It decreases by one at each retraced step.
             */
            unsigned short value{
                static_cast<unsigned short>( count + base )
            };

            /**
             * Retrace the trajectory and fill in cache values for states below `N`.
             */
            for ( unsigned short remaining{ count };
                  remaining != 0;
                  --remaining, --value ) {

                if ( (n.high == 0) && (n.low < N) ) {
                    v[n.low] = value;
                }

                if ( (n.low & std::size_t{ 1 }) == 0 ) {
                    divide_n_by_two( n );
                } else {
                    three_n_plus_one( n );
                }
            }
        }

        /// Print a new record holder.
        if ( curr_max < v[k] ) {
            std::cout << k << "\t" << v[k] << std::endl;
            curr_max = v[k];
        }
    }

    /// Print the largest intermediate value reached.
    if ( max_n.high == 0 ) {
        std::cout << "Maximum n:    0x"
                  << std::hex << max_n.low
                  << std::endl;
    } else {
        std::cout << "Maximum n:    0x"
                  << std::hex << max_n.high
                  << std::setw(16) << std::setfill('0') << max_n.low
                  << std::endl;
    }

    std::cout << " where 2^64 = 0x10000000000000000" << std::endl;

    return 0;
}

/**
 * @brief Divide a two-word integer by two.
 *
 * @param[in,out] n
 * The two-word value to divide.
 *
 * @details
 * If the least significant bit of `n.high` is 1 before the shift,
 * that bit is transferred into the most significant bit of `n.low`.
 */
void divide_n_by_two( double_size_t &n ) {
    bool hi_lsb{ (n.high & 1u) == 1u };

    n.high >>= 1;
    n.low  >>= 1;

    if ( hi_lsb ) {
        n.low |= MSB;
    }
}

/**
 * @brief Compute \(3n+1\) on a two-word unsigned integer.
 *
 * @param[in,out] n
 * The value to update.
 *
 * @details
 * A fast one-word path is used when
 * \f[   3n + 1 \le \texttt{max(std::size_t)}.   \f]
 *
 * Otherwise, explicit two-word arithmetic is used:
 *
 * - first compute \(2n\)
 * - then add \(n\)
 * - then add 1
 */
void three_n_plus_one( double_size_t &n ) {
    constexpr std::size_t max_safe{
        (std::numeric_limits<std::size_t>::max() - 1) / 3
    };

    if ( (n.high == 0) && (n.low <= max_safe) ) {
        n.low = 3*n.low + 1;
    } else {
        std::size_t const low{ n.low };
        std::size_t const high{ n.high };

        std::size_t low2{ low << 1 };
        std::size_t high2{
            (high << 1) | (low >> (8 * sizeof(std::size_t) - 1))
        };

        std::size_t new_low{ low2 + low };
        std::size_t carry{ (new_low < low2) ? std::size_t{ 1 } : 0 };

        std::size_t new_high{ high2 + high + carry };

        ++new_low;

        if ( new_low == 0 ) {
            ++new_high;
        }

        n.high = new_high;
        n.low  = new_low;
    }
}

/**
 * @brief Compute \((3n+1)/2\) on a two-word unsigned integer.
 *
 * @param[in,out] n
 * The value to update.
 *
 * @details
 * For odd \(n\),
 * \f[   \frac{3n+1}{2} = n + \frac{n+1}{2}.   \f]
 *
 * If `n.high == 0` and the final result still fits in one word, the function
 * uses the fast one-word computation
 * \f[   n \leftarrow n + \frac{n+1}{2}.   \f]
 *
 * Otherwise, it computes:
 *
 * - \((n+1)/2\) in two-word arithmetic
 * - adds the original \(n\)
 *
 * @pre
 * This function should be called only when `n` is odd.
 */
void three_n_plus_one_and_divide_by_2( double_size_t &n ) {
    std::size_t const low{ n.low };

    if ( n.high == 0 ) {
        std::size_t const half_up{ (low >> 1) + 1 };

        if ( low <= std::numeric_limits<std::size_t>::max() - half_up ) {
            n.low = low + half_up;
            return;
        }
    }

    std::size_t const high{ n.high };

    /// Compute \((n+1)/2\).
    std::size_t half_low{ (low + 1) >> 1 };
    std::size_t carry{
        (low == std::numeric_limits<std::size_t>::max())
        ? std::size_t{ 1 }
        : 0
    };

    std::size_t half_high{ (high + carry) >> 1 };

    /// Add the original value \(n\).
    std::size_t new_low{ low + half_low };
    std::size_t carry2{ (new_low < low) ? std::size_t{ 1 } : 0 };

    std::size_t new_high{ high + half_high + carry2 };

    n.high = new_high;
    n.low  = new_low;
}
