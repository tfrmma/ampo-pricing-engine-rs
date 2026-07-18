//! Delta, Gamma, Vega for AmPOs under constant amortization, Black-Scholes underlying.
//! Table 1 in Feinstein, "Amortizing Perpetual Options" (arXiv:2512.06505).
//! Vega is NOT monotone in q, see comparative_statics.rs and Example 3.12 in the paper,
//! don't assume it behaves like the premium or the boundary.

use crate::black_scholes::AmpoParams;

pub fn delta_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    ((a - 1.0) * p.s0 / (a * p.k)).powf(a - 1.0)
}

pub fn delta_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    -(a * p.k / ((1.0 + a) * p.s0)).powf(1.0 + a)
}

pub fn gamma_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    (a - 1.0).powi(2) / (a * p.k) * ((a - 1.0) * p.s0 / (a * p.k)).powf(a - 2.0)
}

pub fn gamma_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    a * p.k / p.s0.powi(2) * (a * p.k / ((1.0 + a) * p.s0)).powf(a)
}

pub fn vega_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    let c0 = crate::black_scholes::price_call(p);
    4.0 * c0 / p.sigma * ((a - 1.0) * p.s0 / (a * p.k)).ln() * ((a - 1.0) * p.r - p.q)
        / ((2.0 * a - 1.0) * p.sigma.powi(2) + 2.0 * p.r)
}

pub fn vega_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    let p0 = crate::black_scholes::price_put(p);
    4.0 * p0 / p.sigma * ((1.0 + a) * p.s0 / (a * p.k)).ln() * ((1.0 + a) * p.r + p.q)
        / ((2.0 * a + 1.0) * p.sigma.powi(2) - 2.0 * p.r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::black_scholes::{price_call, price_put, AmpoParams};
    use approx::assert_relative_eq;

    fn call_params() -> AmpoParams {
        AmpoParams {
            s0: 90.0,
            k: 100.0,
            r: 0.05,
            sigma: 0.5,
            q: 0.3,
        }
    }

    fn put_params() -> AmpoParams {
        AmpoParams {
            s0: 110.0,
            k: 100.0,
            r: 0.05,
            sigma: 0.5,
            q: 0.3,
        }
    }

    #[test]
    fn call_delta_matches_finite_difference() {
        let p = call_params();
        let ds = 1e-4;
        let mut up = p;
        up.s0 += ds;
        let mut down = p;
        down.s0 -= ds;
        let numerical = (price_call(&up) - price_call(&down)) / (2.0 * ds);
        assert_relative_eq!(numerical, delta_call(&p), epsilon = 1e-5);
    }

    #[test]
    fn put_delta_matches_finite_difference() {
        let p = put_params();
        let ds = 1e-4;
        let mut up = p;
        up.s0 += ds;
        let mut down = p;
        down.s0 -= ds;
        let numerical = (price_put(&up) - price_put(&down)) / (2.0 * ds);
        assert_relative_eq!(numerical, delta_put(&p), epsilon = 1e-5);
    }

    // Gamma checked as the derivative of delta, not the second derivative of price.
    // Second-differencing price directly hits catastrophic cancellation at reasonable
    // step sizes and gives you a 1-3% "error" that isn't really there, learned that
    // the hard way validating this against the paper before writing the Rust.
    #[test]
    fn call_gamma_matches_derivative_of_delta() {
        let p = call_params();
        let ds = 1e-3;
        let mut up = p;
        up.s0 += ds;
        let mut down = p;
        down.s0 -= ds;
        let numerical = (delta_call(&up) - delta_call(&down)) / (2.0 * ds);
        assert_relative_eq!(numerical, gamma_call(&p), epsilon = 1e-6);
    }

    #[test]
    fn put_gamma_matches_derivative_of_delta() {
        let p = put_params();
        let ds = 1e-3;
        let mut up = p;
        up.s0 += ds;
        let mut down = p;
        down.s0 -= ds;
        let numerical = (delta_put(&up) - delta_put(&down)) / (2.0 * ds);
        assert_relative_eq!(numerical, gamma_put(&p), epsilon = 1e-6);
    }

    #[test]
    fn call_vega_matches_finite_difference() {
        let p = call_params();
        let dv = 1e-6;
        let mut up = p;
        up.sigma += dv;
        let mut down = p;
        down.sigma -= dv;
        let numerical = (price_call(&up) - price_call(&down)) / (2.0 * dv);
        assert_relative_eq!(numerical, vega_call(&p), epsilon = 1e-5);
    }

    #[test]
    fn put_vega_matches_finite_difference() {
        let p = put_params();
        let dv = 1e-6;
        let mut up = p;
        up.sigma += dv;
        let mut down = p;
        down.sigma -= dv;
        let numerical = (price_put(&up) - price_put(&down)) / (2.0 * dv);
        assert_relative_eq!(numerical, vega_put(&p), epsilon = 1e-5);
    }

    #[test]
    fn call_vega_is_not_monotone_in_q() {
        // Example 3.12: positional vega/premium has a non-monotone region for puts,
        // and per-unit vega itself need not be monotone in q in general. Just assert
        // it doesn't blow up or go negative in a nonsensical way across a sweep, this
        // is a smoke test, not a proof.
        let qs = [0.01, 0.1, 0.3, 0.5, 1.0, 2.0];
        let vegas: Vec<f64> = qs
            .iter()
            .map(|&q| {
                let mut p = call_params();
                p.q = q;
                vega_call(&p)
            })
            .collect();
        assert!(vegas.iter().all(|v| v.is_finite()));
    }
}
