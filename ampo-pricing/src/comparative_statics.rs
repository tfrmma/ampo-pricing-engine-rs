//! Sensitivity to the amortization rate q itself, not a hedging Greek. This is for
//! picking canonical q values when setting up a market (Section 4.1: the paper
//! suggests a few fixed rates like 3bps/day, 10bps/day, 50bps/day to segment
//! participants by implied maturity, similar to picking standard expiries).

use crate::black_scholes::AmpoParams;

fn alpha_bar(p: &AmpoParams) -> f64 {
    (p.alpha_call() + p.alpha_put()) / 2.0
}

/// Cor 3.10. Premium is always decreasing in q, both option types.
pub fn dprice_dq_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    let c0 = crate::black_scholes::price_call(p);
    c0 / (p.sigma.powi(2) * alpha_bar(p)) * ((a - 1.0) * p.s0 / (a * p.k)).ln()
}

pub fn dprice_dq_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    let p0 = crate::black_scholes::price_put(p);
    -p0 / (p.sigma.powi(2) * alpha_bar(p)) * ((1.0 + a) * p.s0 / (a * p.k)).ln()
}

/// Cor 3.11. Boundary shrinks toward K as q grows, from above for calls, from below
/// for puts.
pub fn dboundary_dq_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    -p.k / (p.sigma.powi(2) * (a - 1.0).powi(2) * alpha_bar(p))
}

pub fn dboundary_dq_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    p.k / (p.sigma.powi(2) * (1.0 + a).powi(2) * alpha_bar(p))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::black_scholes::{exercise_boundary_call, exercise_boundary_put};
    use approx::assert_relative_eq;

    fn params(q: f64) -> AmpoParams {
        AmpoParams { s0: 90.0, k: 100.0, r: 0.05, sigma: 0.5, q }
    }

    #[test]
    fn call_boundary_derivative_matches_finite_difference() {
        let q0 = 0.5;
        let dq = 1e-6;
        let numerical =
            (exercise_boundary_call(&params(q0 + dq)) - exercise_boundary_call(&params(q0 - dq))) / (2.0 * dq);
        assert_relative_eq!(numerical, dboundary_dq_call(&params(q0)), epsilon = 1e-4);
    }

    #[test]
    fn put_boundary_derivative_matches_finite_difference() {
        let q0 = 0.5;
        let dq = 1e-6;
        let numerical =
            (exercise_boundary_put(&params(q0 + dq)) - exercise_boundary_put(&params(q0 - dq))) / (2.0 * dq);
        assert_relative_eq!(numerical, dboundary_dq_put(&params(q0)), epsilon = 1e-4);
    }

    #[test]
    fn premium_decreasing_in_q_both_types() {
        assert!(dprice_dq_call(&params(0.5)) <= 0.0);
        assert!(dprice_dq_put(&params(0.5)) <= 0.0);
    }

    #[test]
    fn call_boundary_decreasing_toward_strike() {
        assert!(dboundary_dq_call(&params(0.5)) < 0.0);
    }

    #[test]
    fn put_boundary_increasing_toward_strike() {
        assert!(dboundary_dq_put(&params(0.5)) > 0.0);
    }
}
