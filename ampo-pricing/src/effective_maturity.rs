//! Effective maturity: the T for a dated American option (same S0, r, sigma) that
//! matches an AmPO's premium. Purely a pricing equivalence, footnote 3 in the paper
//! is explicit that this does NOT mean the AmPO "expires" at T_eff, it retains most
//! of its notional there. Don't let this leak into anything that treats T_eff as an
//! actual maturity.
//!
//! There's no closed form for a dated American option (calls without dividends never
//! exercise early so European = American there, but puts do), so this pulls in a CRR
//! binomial tree just for this module.

use crate::black_scholes::AmpoParams;

/// CRR binomial tree, American exercise, no dividends. steps=500 is overkill for the
/// precision we need here (this isn't a production pricer, it's a reference point
/// for effective_maturity) but binomial trees are cheap enough that it doesn't matter.
pub fn american_option_price(s0: f64, k: f64, r: f64, sigma: f64, t: f64, is_call: bool, steps: usize) -> f64 {
    if t <= 0.0 {
        return if is_call { (s0 - k).max(0.0) } else { (k - s0).max(0.0) };
    }
    let dt = t / steps as f64;
    let u = (sigma * dt.sqrt()).exp();
    let d = 1.0 / u;
    let disc = (-r * dt).exp();
    let p_up = ((r * dt).exp() - d) / (u - d);

    // terminal payoffs
    let mut values: Vec<f64> = (0..=steps)
        .map(|i| {
            let s_t = s0 * u.powi(i as i32) * d.powi((steps - i) as i32);
            if is_call { (s_t - k).max(0.0) } else { (k - s_t).max(0.0) }
        })
        .collect();

    // backward induction, checking early exercise at every node
    for step in (0..steps).rev() {
        for i in 0..=step {
            let s_node = s0 * u.powi(i as i32) * d.powi((step - i) as i32);
            let continuation = disc * (p_up * values[i + 1] + (1.0 - p_up) * values[i]);
            let intrinsic = if is_call { (s_node - k).max(0.0) } else { (k - s_node).max(0.0) };
            values[i] = continuation.max(intrinsic);
        }
    }
    values[0]
}

/// Bisection for T such that the dated American option price matches the AmPO price.
/// Only defined for q > 0 and t_max large enough that american_option_price(t_max)
/// exceeds the AmPO value, panics otherwise rather than returning a garbage bound.
pub fn effective_maturity_call(p: &AmpoParams, t_max: f64, steps: usize, tol: f64) -> f64 {
    let target = crate::black_scholes::price_call(p);
    bisect_maturity(p.s0, p.k, p.r, p.sigma, true, target, t_max, steps, tol)
}

pub fn effective_maturity_put(p: &AmpoParams, t_max: f64, steps: usize, tol: f64) -> f64 {
    let target = crate::black_scholes::price_put(p);
    bisect_maturity(p.s0, p.k, p.r, p.sigma, false, target, t_max, steps, tol)
}

fn bisect_maturity(
    s0: f64,
    k: f64,
    r: f64,
    sigma: f64,
    is_call: bool,
    target: f64,
    t_max: f64,
    steps: usize,
    tol: f64,
) -> f64 {
    let price_at = |t: f64| american_option_price(s0, k, r, sigma, t, is_call, steps);
    let mut lo = 0.0;
    let mut hi = t_max;
    assert!(
        price_at(hi) >= target,
        "t_max too small, dated option at t_max ({}) doesn't reach the AmPO target ({})",
        price_at(hi),
        target
    );
    while hi - lo > tol {
        let mid = (lo + hi) / 2.0;
        if price_at(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::black_scholes::{price_call, price_put, AmpoParams};
    use approx::assert_relative_eq;

    #[test]
    fn american_call_matches_black_scholes_european() {
        // no dividends, so American call == European call, can cross-check against
        // the standard closed form instead of trusting the tree blindly.
        let (s0, k, r, sigma, t) = (100.0, 100.0, 0.05, 0.3, 1.0);
        let tree_price = american_option_price(s0, k, r, sigma, t, true, 2000);
        let bs_price = black_scholes_european_call(s0, k, r, sigma, t);
        assert_relative_eq!(tree_price, bs_price, epsilon = 5e-3);
    }

    fn black_scholes_european_call(s0: f64, k: f64, r: f64, sigma: f64, t: f64) -> f64 {
        let d1 = ((s0 / k).ln() + (r + 0.5 * sigma.powi(2)) * t) / (sigma * t.sqrt());
        let d2 = d1 - sigma * t.sqrt();
        s0 * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
    }

    // Abramowitz-Stegun approximation, good enough for a cross-check, not for production.
    fn norm_cdf(x: f64) -> f64 {
        0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
    }

    fn erf(x: f64) -> f64 {
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;
        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x = x.abs();
        let t = 1.0 / (1.0 + p * x);
        let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
        sign * y
    }

    #[test]
    fn american_price_increasing_in_maturity() {
        let (s0, k, r, sigma) = (100.0, 100.0, 0.05, 0.3);
        let short = american_option_price(s0, k, r, sigma, 0.5, false, 500);
        let long = american_option_price(s0, k, r, sigma, 2.0, false, 500);
        assert!(long > short);
    }

    #[test]
    fn effective_maturity_recovers_a_price_match_for_call() {
        let p = AmpoParams { s0: 100.0, k: 100.0, r: 0.05, sigma: 0.5, q: 0.3 };
        let ampo_price = price_call(&p);
        let t_eff = effective_maturity_call(&p, 50.0, 500, 1e-3);
        let dated_price = american_option_price(p.s0, p.k, p.r, p.sigma, t_eff, true, 500);
        assert_relative_eq!(dated_price, ampo_price, epsilon = 1e-2);
    }

    #[test]
    fn effective_maturity_recovers_a_price_match_for_put() {
        let p = AmpoParams { s0: 100.0, k: 100.0, r: 0.05, sigma: 0.5, q: 0.3 };
        let ampo_price = price_put(&p);
        let t_eff = effective_maturity_put(&p, 50.0, 500, 1e-3);
        let dated_price = american_option_price(p.s0, p.k, p.r, p.sigma, t_eff, false, 500);
        assert_relative_eq!(dated_price, ampo_price, epsilon = 1e-2);
    }

    // TODO: T_eff should be decreasing in q (higher amortization -> AmPO worth less
    // -> matches a shorter-dated option), haven't added a monotonicity sweep test for
    // that yet, current tests only check a single point matches.
}
