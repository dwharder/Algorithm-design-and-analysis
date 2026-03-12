# Hyperbolic bound for rate-monotonic scheduling

This repository implements a fixed-point integer algorithm for checking the
**Bini–Buttazzo hyperbolic bound** for rate-monotonic scheduling.

The algorithm evaluates the schedulability condition

$$\prod_{k=1}^{n} (1 + u_k) \le 2$$

where $u_k$ is the utilization of task $k$.

The implementation is designed for **embedded systems** where:

- floating-point arithmetic may not be available
- 64-bit arithmetic may be expensive
- deterministic integer behavior is required

All calculations are performed using **32-bit integers only**. Each step of the computation rounds **upward**, ensuring that the stored product is always an **upper bound** on the true value of

$$\prod_{k=1}^{n} (1 + u_k).$$

The maintained invariant is that the stored value is always **greater than or equal to the true product**.

Consequently, the implementation is **conservative**: if the algorithm accepts a task with a given utilization factor, then the hyperbolic schedulability condition above is guaranteed to be satisfied.

This guarantee cannot be obtained using standard floating-point arithmetic, which typically rounds to the nearest representable value, nor with naive integer arithmetic that truncates fractional values. Both approaches may round downward and therefore risk producing an underestimate of the true product.

The hyperbolic bound itself follows from the result of **Bini and Buttazzo** that a set of periodic tasks scheduled by rate-monotonic priority assignment is schedulable if

$$ \prod_{k=1}^n (1 + u_k) \le 2,$$

where $u_k$ is the utilization of task $k$. This condition is sufficient (though not necessary) and can be evaluated incrementally as tasks are added.

---

# Representation

Each factor $1 + u_k$ is stored using a **32-bit binary mantissa** encoding a value of the form

$$1 + \frac{\text{mantissa}}{2^{32}},$$

chosen so that the encoded value is **greater than or equal to** the true factor $1 + u_k$. Thus each stored factor is an upward-rounded approximation of the true value.
