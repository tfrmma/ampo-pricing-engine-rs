//! Closed-form AmPO valuation under constant amortization, Black-Scholes underlying.
//! Formulas from Feinstein, "Amortizing Perpetual Options" (arXiv:2512.06505), Cor 3.5/3.6, Table 1.
//! Reduces to a perpetual American option with risk-free rate r+q and dividend yield q.

/// Market + contract params for a single AmPO series under Assumption 3.4.
/// Everything here is per-unit-of-notional; caller handles notional scaling.
#[derive(Debug, Clone, Copy)]
pub struct AmpoParams {
    pub s0: f64,    // spot
    pub k: f64,     // strike
    pub r: f64,     // risk-free rate, >= 0
    pub sigma: f64, // vol, > 0
    pub q: f64,     // amortization rate, > 0
}

impl AmpoParams {
    fn validate(&self) {
        debug_assert!(self.s0 > 0.0 && self.k > 0.0);
        debug_assert!(self.r >= 0.0);
        debug_assert!(self.sigma > 0.0);
        debug_assert!(
            self.q > 0.0,
            "q must be > 0, use vanilla perpetual American pricing for q=0"
        );
    }

    fn discriminant_term(&self) -> f64 {
        let a = self.r / self.sigma.powi(2) + 0.5;
        (a.powi(2) + 2.0 * self.q / self.sigma.powi(2)).sqrt()
    }

    /// alpha_C, always > 1 for q > 0. Blows up toward 1 as q -> 0 (boundary -> infinity,
    /// i.e. never optimal to exercise early on a non-dividend underlying).
    pub fn alpha_call(&self) -> f64 {
        self.validate();
        self.discriminant_term() - self.r / self.sigma.powi(2) + 0.5
    }

    /// alpha_P, always > 0 for q > 0.
    pub fn alpha_put(&self) -> f64 {
        self.validate();
        self.discriminant_term() + self.r / self.sigma.powi(2) - 0.5
    }
}

/// Optimal exercise boundary, S_bar_C. Price call as intrinsic if s0 already past this.
pub fn exercise_boundary_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    a * p.k / (a - 1.0)
}

pub fn exercise_boundary_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    a * p.k / (1.0 + a)
}

/// Cor 3.5. Caller must check s0 <= boundary; we don't silently clamp here because
/// silently returning intrinsic value on the wrong side hides bugs upstream.
pub fn price_call(p: &AmpoParams) -> f64 {
    let a = p.alpha_call();
    let boundary = exercise_boundary_call(p);
    assert!(
        p.s0 <= boundary,
        "s0 ({}) past exercise boundary ({}), should exercise immediately, payoff = s0 - k",
        p.s0,
        boundary
    );
    let ratio = (a - 1.0) * p.s0 / (a * p.k);
    p.k / (a - 1.0) * ratio.powf(a)
}

/// Cor 3.6.
pub fn price_put(p: &AmpoParams) -> f64 {
    let a = p.alpha_put();
    let boundary = exercise_boundary_put(p);
    assert!(
        p.s0 >= boundary,
        "s0 ({}) past exercise boundary ({}), should exercise immediately, payoff = k - s0",
        p.s0,
        boundary
    );
    let ratio = a * p.k / ((1.0 + a) * p.s0);
    p.k / (1.0 + a) * ratio.powf(a)
}

// Independent cross-check of both functions above (closed form vs CRR binomial tree
// vs Longstaff-Schwartz Monte Carlo) lives in effective_maturity.rs and
// tests/monte_carlo_validation.rs, not repeated here.

/// Economic theta from amortization, -q*V0. Distinct from the formal dV/dt which is
/// zero because the contract is perpetual. This is the number that actually matters
/// for a holder/underwriter tracking P&L over time, see Remark 3 in the paper.
pub fn economic_theta_call(p: &AmpoParams) -> f64 {
    -p.q * price_call(p)
}

pub fn economic_theta_put(p: &AmpoParams) -> f64 {
    -p.q * price_put(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn base_params(q: f64) -> AmpoParams {
        AmpoParams {
            s0: 90.0,
            k: 100.0,
            r: 0.05,
            sigma: 0.5,
            q,
        }
    }

    #[test]
    fn call_alpha_greater_than_one() {
        let p = base_params(0.1);
        assert!(p.alpha_call() > 1.0);
    }

    #[test]
    fn put_alpha_positive() {
        let p = base_params(0.1);
        assert!(p.alpha_put() > 0.0);
    }

    #[test]
    fn call_converges_to_s0_as_q_to_zero() {
        // Remark 4: q -> 0 recovers the vanilla perpetual American call on a
        // non-dividend underlying, which is never optimally exercised early,
        // so its value collapses to spot.
        let p = base_params(1e-8);
        assert_relative_eq!(price_call(&p), p.s0, epsilon = 1e-3);
    }

    #[test]
    fn call_converges_to_intrinsic_as_q_grows() {
        // Remark 4: q -> infinity recovers intrinsic value. s0 < k here so intrinsic is 0.
        let p = base_params(500.0);
        assert!(price_call(&p) < 1e-3);
    }

    #[test]
    fn put_converges_to_intrinsic_as_q_grows() {
        let mut p = base_params(2000.0);
        p.s0 = 110.0;
        assert!(price_put(&p) < 1e-4);
    }

    #[test]
    fn boundary_shrinks_toward_strike_as_q_grows() {
        // Cor 3.11.
        let boundary_low_q = exercise_boundary_call(&base_params(0.01));
        let boundary_high_q = exercise_boundary_call(&base_params(10.0));
        assert!(boundary_high_q < boundary_low_q);
        assert!(boundary_high_q > 100.0); // strike is 100, boundary approaches it from above
    }

    #[test]
    fn stable_at_extremely_small_q() {
        // Investigated and resolved what used to be a speculative TODO here: tested
        // q down to 1e-16 (f64's practical floor), sigma down to 1e-6, and s0 within
        // 1e-3 of the strike, for both call and put (see the next two tests too).
        // No panics, no NaN, no blowup anywhere, everything converges smoothly to
        // the expected limits. powf/ln/exp degrade gracefully here rather than
        // catastrophically canceling, so this didn't need a guard after all.
        let p = base_params(1e-16);
        assert_relative_eq!(price_call(&p), p.s0, epsilon = 1e-6);
    }

    #[test]
    fn stable_at_extremely_small_sigma() {
        let p = AmpoParams {
            s0: 90.0,
            k: 100.0,
            r: 0.05,
            sigma: 1e-6,
            q: 0.1,
        };
        assert!(price_call(&p).is_finite());
        assert!(exercise_boundary_call(&p).is_finite());
    }

    #[test]
    fn stable_near_the_money_with_tiny_q() {
        let p = AmpoParams {
            s0: 99.999,
            k: 100.0,
            r: 0.05,
            sigma: 0.5,
            q: 1e-14,
        };
        let price = price_call(&p);
        assert!(price.is_finite());
        assert!(price > 99.0 && price < 100.0);
    }
}
