#ifndef CA_UWATERLOO_DWHARDER_HYPERBOLIC_BOUND
#define CA_UWATERLOO_DWHARDER_HYPERBOLIC_BOUND

#include <cstdint>

/**
 * @file Hyperbolic_bound.cpp
 *
 * @brief Integer implementation of the Bini–Buttazzo hyperbolic schedulability bound.
 *
 * @author Douglas Wilhelm Harder
 *
 * @details
 * This file implements the schedulability test introduced by
 * Bini and Buttazzo for rate-monotonic scheduling.
 *
 * A task set with utilizations \(u_k = C_k/T_k\) satisfies the bound
 *
 * \f[   \prod_{k = 1}^N (1 + u_k) \le 2 .   \f]
 *
 * If the bound holds, the task set is guaranteed to be schedulable
 * under rate-monotonic priorities.
 *
 * This implementation avoids floating-point arithmetic entirely and
 * performs all calculations using 32-bit integers.  The product
 *
 * \f[   \prod_{k = 1}^N (1 + u_k)   \f]
 *
 * is represented using a 32-bit binary mantissa.  All rounding is
 * performed upward so that the stored value is always a conservative
 * over-approximation of the true product.  Consequently, the algorithm
 * never incorrectly declares an unschedulable system schedulable.
 */


/**
 * @brief Maintains the running product used in the Bini–Buttazzo hyperbolic bound.
 *
 * @details
 * The stored value represents
 *
 * \f[   1 + \text{mantissa}/2^{32}   \f]
 *
 * and encodes the running product
 *
 * \f[   \prod_{k = 1}^N (1 + u_k).   \f]
 *
 * The internal encoding stores only the fractional mantissa.
 *
 * The special encoding
 *
 * \f[   \texttt{mantissa\_} = 0   \f]
 *
 * represents the exact value \(2.0\).
 *
 * This value indicates that the schedulability bound has been
 * reached and that no further tasks can be added.
 *
 * All arithmetic is performed so that
 *
 * \f[   \text{true value} \le \text{stored value}.   \f]
 *
 * Thus rounding is always conservative.
 *
 * This class assumes that at least one base periodic task exists,
 * so the initial product is always strictly greater than \(1\).
 */

class Hyperbolic_bound {
    static_assert( sizeof( std::uint32_t ) == 4,
                  "Hyperbolic_bound requires 32-bit integers" );

    using u32 = std::uint32_t;

    public:

        /**
         * @brief Construct the bound using the base task utilization factor.
         *
         * @param base_task_utilization_factor
         * Encoded utilization factor for the base periodic task.
         *
         * @details
         * The initial value represents
         *
         * \f[   1 + u_{\text{base}} .   \f]
         *
         * If the supplied value is zero, it is replaced with the smallest
         * non-zero encoding so that the stored product always represents
         * a value strictly greater than \(1\).
         */

        explicit Hyperbolic_bound( u32 base_task_utilization_factor )
        : mantissa_{ base_task_utilization_factor ? base_task_utilization_factor : 1 } {
            // Initialization complete
        }

        /**
         * @brief Attempt to add a new periodic task.
         *
         * @param utilization_factor
         * Encoded utilization factor \(1 + u_k\).
         *
         * @return true if the task can be added while maintaining
         *         \f$ \prod (1 + u_k) \le 2 \f$.
         *
         * @return false if adding the task would cause the product
         *         to exceed \(2\).  In this case the stored value is
         *         left unchanged.
         *
         * @details
         * This routine multiplies
         *
         * \f[   (1 + m)(1 + u)   \f]
         *
         * using only 32-bit arithmetic.
         *
         * The mantissas are split into 16-bit halves so that the
         * multiplication can be performed using four 16×16 products.
         *
         * The discarded lower bits are examined and the result is
         * rounded upward whenever necessary so that the stored value
         * remains a conservative over-approximation.
         */

        [[nodiscard]] bool add_task( u32 utilization_factor ) noexcept {

            // If either the current product is 2 or the utlization factor is 1,
            //    we cannot add this next task.
            if ( (mantissa_ == 0) || (utilization_factor == 0) ) {
                return false;
            }

            //       1.m0 m1
            //     x 1.u0 u1
            //       -------------
            //               u1*m1
            //            u1*m0
            //            u1
            //            u0*m1
            //         u0*m0
            //         u0
            //     + 1.m0 m1
            //       -------------
            //                 =p23==
            //             =p12a=
            //             u1
            //             =p12b=
            //         =p01==
            //         u0
            //     + 1.m0  m1
            //       -------------
            //                 p2c s3
            //             p1c p2b
            //             u1
            //             p1b p2a
            //         p0  p1a
            //         u0
            //     + 1.m0  m1
            //       -------------
            //       1.s0  s1  s2  s3
            //
            // \note
            // We could reduce the four products to three using a
            // Karatsuba-style identity, but the intermediate product
            // (m0 + m1)*(u0 + u1) may not fit in 32 bits.

            u32 m1{ mantissa_ & 0xffffu };
            u32 m0{ mantissa_ >> 16 };

            u32 u1{ utilization_factor & 0xffffu };
            u32 u0{ utilization_factor >> 16 };

            u32 p01 { m0*u0 };
            u32 p12a{ m0*u1 };
            u32 p12b{ m1*u0 };
            u32 p23 { m1*u1 };

            u32 p0  { p01 >> 16 };
            u32 p1a { p01 & 0xffffu };
            u32 p1b { p12a >> 16 };
            u32 p1c { p12b >> 16 };
            u32 p2a { p12a & 0xffffu };
            u32 p2b { p12b & 0xffffu };
            u32 p2c { p23 >> 16 };

            u32 s3  { p23 & 0xffffu };

            u32 s2{ p2a + p2b + p2c };

            u32 s1{ p1a + p1b + p1c + m1 + u1 + (s2 >> 16) };
            s2 &= 0xffffu;

            if ( (s3 > 0) || (s2 > 0) ) {
                ++s1;
            }

            u32 s0{ p0 + m0 + u0 + (s1 >> 16) };
            s1 &= 0xffffu;

            if ( (s0 > 0x10000u) || ((s0 == 0x10000u) && (s1 != 0)) ) {
                return false;
            }

            /// \invariant The stored value is always an upper bound on the true product.
            mantissa_ = (s0 << 16) | s1;
            return true;
        }

        /**
         * @brief Remove an existing periodic task.
         *
         * @param utilization_factor
         * Encoded utilization factor \(1 + u_k\).
         *
         * @return true if the task can be added while maintaining
         *         \f$ \prod (1 + u_k) \le 2 \f$.
         *
         * @return false if adding the task would cause the product
         *         to exceed \(2\).  In this case the stored value is
         *         left unchanged.
         *
         * @details
         * This routine multiplies
         *
         * \f[   (1 + m)(1 + u)   \f]
         *
         * using only 32-bit arithmetic.
         *
         * The mantissas are split into 16-bit halves so that the
         * multiplication can be performed using four 16×16 products.
         *
         * The discarded lower bits are examined and the result is
         * rounded upward whenever necessary so that the stored value
         * remains a conservative over-approximation.
         *
         * \warning
         * This operation should be used only occasionally. Repeated sequences
         * of adding and removing tasks will accumulate rounding error.
         * If tasks are frequently inserted or removed, it is better to
         * recompute the bound from the full set of tasks.
         */

        void remove_task( u32 utilization_factor ) noexcept {
            // If the utilization is 1.00000000, set the mantissa
            // to be 1.00000001, the smallest product other than 1.
            if ( utilization_factor == 0 ) {
                mantissa_ = 1;
                return;
            }

            //
            // Numerator N is the current value:
            //
            //   N = 1.mantissa_
            //
            // except that mantissa_ == 0 represents exactly 2.0.
            //
            // Denominator D is:
            //
            //   D = 1.utilization_factor
            //
            // We want
            //
            //   Q = ceil( N / D )
            //
            // encoded as a mantissa.  Since valid removal should produce
            // a value in (1, 2), we write
            //
            //   Q = 1 + (N - D)/D.
            //
            // Thus we first subtract D once, and then compute the 32-bit
            // fractional expansion of (N - D)/D, rounding upward if any
            // remainder remains.
            //

            u32 numerator_hi{ 1 };
            u32 numerator_lo{ mantissa_ };

            if ( mantissa_ == 0 ) {
                // exactly 2.0
                numerator_hi = 2;
                numerator_lo = 0;
            }

            u32 denom_hi{ 1 };
            u32 denom_lo{ utilization_factor };

            //
            // If N <= D, then removing this factor would produce a value
            // less than or equal to 1.0.  In this design, clamp to the
            // smallest representable value greater than 1.
            //
            if (
                (numerator_hi < denom_hi)
                || ((numerator_hi == denom_hi) && (numerator_lo <= denom_lo))
            ) {
                mantissa_ = 1;
                return;
            }

            //
            // Compute R = N - D.
            //
            {
                u32 old_lo{ numerator_lo };
                numerator_lo -= denom_lo;
                numerator_hi -= denom_hi + (old_lo < denom_lo ? 1u : 0);
            }

            //
            // Now compute the mantissa of R / D.
            //
            mantissa_ = 0;

            for ( u32 mask{ 0x80000000u }; mask; mask >>= 1 ) {
                // Shift remainder left by one bit.
                numerator_hi = (numerator_hi << 1) | (numerator_lo >> 31);
                numerator_lo <<= 1;

                // If remainder >= denominator, subtract and set this bit.
                if (
                    (numerator_hi > denom_hi)
                    || ((numerator_hi == denom_hi) && (numerator_lo >= denom_lo))
                ) {

                    u32 old_lo{ numerator_lo };
                    numerator_lo -= denom_lo;
                    numerator_hi -= denom_hi + (old_lo < denom_lo ? 1u : 0);

                    mantissa_ |= mask;
                }
            }

            // Round upward if a remainder remains.
            if ( (numerator_hi != 0) || (numerator_lo != 0) ) {
                ++mantissa_;

                //
                // If this overflows from 0xffffffff to 0, that correctly
                // represents exactly 2.0 in this encoding.
                //
            }

            //
            // Preserve the convention that the value is strictly greater than 1.
            // If the result is exactly 1.0, represent it by the smallest
            // positive mantissa instead.
            //
            if ( mantissa_ == 0x00000000u ) {
                // Leave as 0 only if this came from upward rounding to 2.0.
                // Otherwise, if you never want exact 2.0 here, adjust policy.
            }
        }

        /**
         * @brief Return the raw 32-bit encoding of the product.
         *
         * @return 32-bit mantissa encoding the current product.
         */

        [[nodiscard]] constexpr u32 raw_encoding() const noexcept {
            return mantissa_;
        }

        /**
         * @brief Reset the bound to the utilization factor of a base task.
         *
         * @param base_task_utilization_factor
         * Encoded utilization factor representing
         * \f$1 + u_{\text{base}}\f$.
         *
         * @details
         * This function reinitializes the running product used in the
         * Bini–Buttazzo hyperbolic schedulability bound:
         *
         * \f[   P = \prod_{k = 1}^N (1 + u_k).   \f]
         *
         * The new value of the product becomes
         *
         * \f[   P = 1 + u_{\text{base}}.   \f]
         *
         * A value of zero is treated as invalid input and replaced with the
         * smallest non-zero encoding so that the stored value always represents
         * a number strictly greater than \f$1\f$.
         *
         * This class assumes that the system always contains at least one
         * periodic task, so the value \f$1\f$ is never represented explicitly.
         */

        void reset( u32 base_task_utilization_factor ) noexcept {
            mantissa_ = base_task_utilization_factor ? base_task_utilization_factor : 1;
        }

        /**
         * @brief Recompute the hyperbolic bound from a range of task factors.
         *
         * @param first Pointer to the first encoded utilization factor.
         * @param last  Pointer one past the final encoded utilization factor.
         *
         * @return
         * Pointer to the first factor that would cause the product to exceed
         * \f$2\f$, or `last` if all factors are accepted.
         *
         * @details
         * This routine recomputes
         *
         * \f[   P = \prod_{k = 1}^N (1 + u_k)   \f]
         *
         * from the sequence of encoded utilization factors stored in the range
         * \f$[\,\texttt{first},\,\texttt{last})\f$.
         *
         * The first factor initializes the product, and subsequent factors are
         * incorporated using `add_task()`.
         *
         * If incorporating a factor would cause
         *
         * \f[   P > 2 ,   \f]
         *
         * the stored value remains unchanged and a pointer to the offending
         * factor is returned.
         *
         * If the entire range satisfies the bound, the function returns `last`.
         *
         * If the range is empty, the stored value is reset to the smallest
         * representable value greater than \f$1\f$.
         */

        [[nodiscard]] u32 const *recompute( u32 const *first, u32 const *last ) noexcept {
            if ( first == last ) {
                mantissa_ = 1;
                return last;
            }

            mantissa_ = (*first != 0) ? *first : 1;
            ++first;

            for ( ; first != last; ++first ) {
                if ( !add_task( *first ) ) {
                    return first;
                }
            }

            return last;
        }

        /**
         * @brief Encode a task utilization as a 32-bit binary fraction.
         *
         * @param usage  Worst-case execution time.
         * @param period Task period.
         *
         * @return Encoded mantissa representing
         *
         * \f[   1 + \left\lceil 2^{32} \frac{\text{usage}}{\text{period}} \right\rceil 2^{-32}.   \f]
         *
         * @details
         * This routine computes the binary expansion of
         *
         * \f[   \frac{\text{usage}}{\text{period}}   \f]
         *
         * using integer arithmetic only.
         *
         * The expansion is truncated after 32 bits and rounded upward
         * whenever a remainder remains.
         *
         * Invalid inputs (\f$ usage = 0 \f$ or \f$ usage \ge period \f$)
         * return zero.
         */

        static constexpr u32 encode_utilization_factor( u32 usage, u32 period ) noexcept {
            if ( (usage == 0) || (usage >= period) ) {
                return 0;
            }

            u32 mantissa{ 0 };

            for ( u32 mask{ 0x80000000u }; mask && (usage != 0); mask >>= 1 ) {
                usage <<= 1;

                if ( usage >= period ) {
                    mantissa |= mask;
                    usage -= period;
                }
            }

            if ( usage > 0 ) {
                ++mantissa;
            }

            return mantissa;
        }

    private:

        /**
         * @brief 32-bit mantissa encoding the current product.
         *
         * @details
         * The stored value represents
         *
         * \f[   1 + \text{mantissa\_}/2^{32}.   \f]
         *
         * The value `0` represents exactly \(2.0\).
         */

        u32 mantissa_;
};

#endif
