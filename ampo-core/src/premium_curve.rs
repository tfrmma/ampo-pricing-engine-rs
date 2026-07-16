//! Premium curves as functions of utilization U = X/C. Definition 3.1 requires
//! non-decreasing, positive on (0,1], with a specific limiting behavior as U -> 1:
//! infinity for calls (can't fully collateralize unbounded upside), K for puts
//! (bounded by the strike). Definition 3.3 pairs each premium function P with its
//! net premium phi = integral of P, Prop 3.4 says this pairing is a bijection under
//! mild regularity. The two curves here are the concrete examples from the paper
//! (Example 3.2/3.5), not the only valid choices, just the ones with a documented
//! closed form.

pub trait PremiumFunction {
    /// P(U), marginal cost of the "next" unit of notional.
    fn premium(&self, u: f64) -> f64;
    /// phi(U) = integral_0^U P(u) du, realized cost of the position.
    fn net_premium(&self, u: f64) -> f64;
}

pub struct CallPremiumCurve;

impl PremiumFunction for CallPremiumCurve {
    fn premium(&self, u: f64) -> f64 {
        debug_assert!((0.0..=1.0).contains(&u), "call utilization must be in [0,1]");
        if u >= 1.0 {
            return f64::INFINITY;
        }
        2.0 * u / (1.0 - u).powi(3)
    }

    fn net_premium(&self, u: f64) -> f64 {
        debug_assert!((0.0..=1.0).contains(&u));
        if u >= 1.0 {
            return f64::INFINITY;
        }
        (u / (1.0 - u)).powi(2)
    }
}

pub struct PutPremiumCurve {
    pub k: f64,
}

impl PremiumFunction for PutPremiumCurve {
    fn premium(&self, u: f64) -> f64 {
        debug_assert!((0.0..=1.0).contains(&u), "put utilization must be in [0,1]");
        self.k * u
    }

    fn net_premium(&self, u: f64) -> f64 {
        debug_assert!((0.0..=1.0).contains(&u));
        self.k * u.powi(2) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Prop 3.4: phi'(U) == P(U). Check by finite difference instead of trusting the
    // closed forms, this is exactly the kind of thing where a sign or exponent typo
    // is easy to make and easy to miss by eye.
    #[test]
    fn call_net_premium_derivative_matches_premium() {
        let curve = CallPremiumCurve;
        for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let du = 1e-6;
            let numerical = (curve.net_premium(u + du) - curve.net_premium(u - du)) / (2.0 * du);
            assert_relative_eq!(numerical, curve.premium(u), epsilon = 1e-4);
        }
    }

    #[test]
    fn put_net_premium_derivative_matches_premium() {
        let curve = PutPremiumCurve { k: 100.0 };
        for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let du = 1e-6;
            let numerical = (curve.net_premium(u + du) - curve.net_premium(u - du)) / (2.0 * du);
            assert_relative_eq!(numerical, curve.premium(u), epsilon = 1e-4);
        }
    }

    #[test]
    fn call_net_premium_blows_up_toward_full_utilization() {
        let curve = CallPremiumCurve;
        assert!(curve.net_premium(0.999) > curve.net_premium(0.9));
        assert!(curve.net_premium(0.9999) > 1000.0);
    }

    #[test]
    fn put_net_premium_bounded_by_half_strike() {
        // phi_put(1) = K/2, this is the put's finite limit, contrast with the call
        // which diverges. Both are consistent with Def 3.1's boundary conditions on
        // P (call diverges, put -> K), integrated once more.
        let curve = PutPremiumCurve { k: 100.0 };
        assert_relative_eq!(curve.net_premium(1.0), 50.0);
    }

    #[test]
    fn call_premium_is_infinite_at_full_utilization() {
        // regression test: this used to panic instead of returning the limit
        // Def 3.1 requires, caught by the invariants fuzz test hitting u==1.0
        // exactly via a withdrawal that fully drained collateral headroom.
        let curve = CallPremiumCurve;
        assert!(curve.premium(1.0).is_infinite());
        assert!(curve.net_premium(1.0).is_infinite());
    }

    #[test]
    fn call_premium_zero_at_zero_utilization() {
        let curve = CallPremiumCurve;
        assert_relative_eq!(curve.net_premium(0.0), 0.0);
    }
}
