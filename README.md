# ampo-pricing-engine

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
diagnosed ~3-7% low bias specifically on calls. Investigated across three
independent axes before accepting this as a structural limitation rather than a
bug: time discretization (200→2000 steps, unchanged), regression basis
(degree 2→5, unchanged, no trend), and simulation horizon (T=10→80y with dt held
fixed, plateaus around 3-4%, doesn't converge to zero). A parallel isolated check
(European-only pricing, no exercise decision) matches its independently-computed
reference tightly, ruling out a path-simulation bug. What's left is the standard
regression-based LSM exercise policy being a well-documented biased-low estimator,
worst when almost all of an option's value comes from near-immediate exercise
(which is exactly this regime: large q makes early exercise dominant). Closing
this for real needs a different method (e.g. Andersen-Broadie primal-dual bounds),
which is a bigger undertaking than a cross-check module warrants. The CRR binomial
tree in `effective_maturity.rs` doesn't have this problem (exact backward
induction, no regression) and is the trustworthy independent check for the call
boundary. Tolerances are set from the measured bias, not picked to force a pass.
Full diagnostic trail reproducible via `cargo test -- --ignored --nocapture`.

### Numerical stability at extreme parameters: investigated, not an issue

Initially flagged as a TODO (denominator `alpha_C - 1 -> 0` as `q -> 0` looked like
a plausible catastrophic-cancellation risk). Tested directly: `q` down to `1e-16`
(f64's practical floor), `sigma` down to `1e-6`, `s0` within `1e-3` of the strike,
both option types. No panics, no NaN, no blowup anywhere, everything converges
smoothly to the expected limits. `powf`/`ln`/`exp` evidently degrade gracefully
here. Kept as regression tests in `black_scholes.rs` rather than as a lingering
unresolved TODO.

## Running

```
cargo test --workspace --release
```

63 unit tests, 3 integration tests, 1 ignored diagnostic test (run with
`cargo test -- --ignored --nocapture` to see the bias diagnosis reproduced).
The Monte Carlo tests are meaningfully slower than the rest (~5-10s), hence
`--release`.
