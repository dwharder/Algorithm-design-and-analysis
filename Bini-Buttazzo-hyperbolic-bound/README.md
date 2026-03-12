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

All calculations are performed using **32-bit integers only**.

---

# Representation

Each factor $1 + u_k$ is stored using a **32-bit binary mantissa**
