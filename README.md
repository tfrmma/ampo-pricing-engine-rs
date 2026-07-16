# ampo-pricing-engine-rs

Reference implementation of **Amortizing Perpetual Options (AmPOs)**: options with
exponentially decaying notional and no fixed expiry, priced in closed form and
traded in an on-chain market with formally provable solvency and path-independence
invariants.

Based on two papers:

- Feinstein, Z. (2026). *Amortizing Perpetual Options*. [arXiv:2512.06505](https://arxiv.org/abs/2512.06505)
  — the instrument and its closed-form valuation under Black-Scholes.
- Bichuch, M. & Feinstein, Z. (2026). *Designing On-Chain Options*. [arXiv:2605.19146](https://arxiv.org/abs/2605.19146)
  — the market mechanism (utilization-based premium curves, the six market
  operations, solvency/path-independence proofs) and two applications.

Kept separate from [`options-pricing-engine-rs`](https://github.com/tfrmma/options-pricing-engine-rs):
different problem family (reduced-form, aggregate-risk mutualization vs. single-asset
stochastic vol pricing), different identity. `ampo-core` and `ampo-applications` have
no dependency on that repo at all; a future `ampo-pricing` Greeks refactor may pull in
its dual-number AD (`ad.rs`) instead of hand-rolled closed forms, not done yet.

## Layout

```
ampo-pricing/        the instrument's fair value (2512.06505)
  black_scholes.rs      closed-form price, exercise boundary (Cor 3.5, 3.6)
  greeks.rs              Delta, Gamma, Vega (Table 1)
  comparative_statics.rs sensitivity to q itself (Cor 3.10, 3.11)
  effective_maturity.rs  CRR binomial tree + bisection, maps an AmPO to the
                         dated American option it's priced equivalent to
  tests/monte_carlo_validation.rs
                         independent Longstaff-Schwartz cross-check, see below

ampo-core/            the market mechanism (2605.19146, Section 3)
  payoff.rs              contract definition, notional decay, physical settlement
  premium_curve.rs       P(U), phi(U), the two worked examples from the paper
  market.rs              total premium function Phi(X,C)
  operations.rs          the six market operations, all as delta_phi
  invariants.rs          fuzz-style tests of solvency and path independence

ampo-applications/    two worked applications (2605.19146, Section 4)
  protective_put.rs      endogenous, oracle-free lending collateral
  depeg_insurance.rs     PSM-as-explicit-option
```

## Status

`ampo-pricing` and `ampo-core` are complete against what the two papers specify.
`ampo-applications` covers the two applications the design paper describes; it's a
thin wrapper reusing `ampo-core`, not new math.

Every closed-form formula here was independently checked against finite differences
(and, for pricing, against a CRR binomial tree and a Longstaff-Schwartz Monte Carlo)
before being trusted. Two bugs were caught this way during development, not by
inspection: a numerical instability at full utilization in `CallPremiumCurve` (fixed,
regression test in `premium_curve.rs`), and a test-tolerance issue mistaken for a
pricing bug in early `black_scholes.rs` development (fixed, was the test, not the
formula, verified independently in Python first).

### Known limitation: LSM low bias on calls

The Longstaff-Schwartz cross-check in `monte_carlo_validation.rs` has a genuine,
diagnosed ~3-7% low bias specifically on calls, confirmed to be neither a
discretization artifact (doesn't shrink from 200 to 2000 time steps) nor a
path-simulation bug (an isolated European-only simulation matches its closed-form
reference tightly). It's the standard regression-based LSM exercise policy being
measurably suboptimal when almost all of an option's value comes from near-immediate
exercise, which is exactly the regime these examples are in (effective maturity
~1.3y against a 10y simulation horizon). The CRR binomial tree in
`effective_maturity.rs` doesn't have this problem (exact backward induction, no
regression) and is the trustworthy independent check for the call boundary. The MC
test tolerances are set from the empirically measured bias, not picked to force a
pass, see the module doc comment for the full diagnostic trail.

### Toolchain note

This was built against rustc/cargo 1.75 (Ubuntu apt default). `proptest`'s current
release pulls a `getrandom` version requiring `edition2024` (cargo 1.85+), which
doesn't resolve here, so `ampo-core/invariants.rs` uses a hand-rolled xorshift64 PRNG
instead of `proptest`. Functionally equivalent for this use case, but doesn't get
proptest's shrinking on failure. Worth revisiting with a current toolchain.

## Running

```
cargo test --workspace --release
```

57 unit tests, 3 integration tests, 1 ignored diagnostic test (run with
`cargo test -- --ignored --nocapture` to see the bias diagnosis reproduced).
The Monte Carlo tests are meaningfully slower than the rest (~5-10s), hence
`--release`.
