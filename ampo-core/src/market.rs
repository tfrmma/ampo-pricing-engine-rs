//! Total premium function Phi(X,C) and market state. Def 3.6: Phi(X,C) = C*phi(X/C)
//! for X <= C, infinity otherwise, with Phi(0,0) := 0 by continuity. Prop 3.7 proves
//! Phi is strictly increasing in X, non-increasing in C, convex, positive homogeneous.
//! We don't re-derive those properties here, invariants.rs tests the consequences
//! that actually matter operationally (solvency, path independence).

use crate::premium_curve::PremiumFunction;

/// Open interest notional and collateral for one AmPO series. Both must stay
/// non-negative and X <= C is the solvency constraint, Phi returns infinity if
/// violated rather than a garbage number.
#[derive(Debug, Clone, Copy)]
pub struct MarketState {
    pub x: f64,
    pub c: f64,
}

pub fn total_premium(curve: &impl PremiumFunction, state: MarketState) -> f64 {
    let MarketState { x, c } = state;
    debug_assert!(x >= 0.0 && c >= 0.0);
    if x > c {
        return f64::INFINITY;
    }
    if c == 0.0 {
        return 0.0; // Phi(0,0) := 0
    }
    c * curve.net_premium(x / c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::premium_curve::{CallPremiumCurve, PutPremiumCurve};
    use approx::assert_relative_eq;

    #[test]
    fn phi_zero_at_origin() {
        let curve = CallPremiumCurve;
        assert_relative_eq!(total_premium(&curve, MarketState { x: 0.0, c: 0.0 }), 0.0);
    }

    #[test]
    fn phi_infinite_when_over_utilized() {
        let curve = CallPremiumCurve;
        assert!(total_premium(&curve, MarketState { x: 10.0, c: 5.0 }).is_infinite());
    }

    #[test]
    fn phi_zero_when_no_open_interest() {
        let curve = PutPremiumCurve { k: 100.0 };
        assert_relative_eq!(total_premium(&curve, MarketState { x: 0.0, c: 50.0 }), 0.0);
    }

    // Prop 3.7: strictly increasing in X, holding C fixed.
    #[test]
    fn phi_strictly_increasing_in_x() {
        let curve = CallPremiumCurve;
        let c = 100.0;
        let phi_lo = total_premium(&curve, MarketState { x: 20.0, c });
        let phi_hi = total_premium(&curve, MarketState { x: 40.0, c });
        assert!(phi_hi > phi_lo);
    }

    // Prop 3.7: non-increasing in C, holding X fixed.
    #[test]
    fn phi_non_increasing_in_c() {
        let curve = CallPremiumCurve;
        let x = 20.0;
        let phi_lo_c = total_premium(&curve, MarketState { x, c: 40.0 });
        let phi_hi_c = total_premium(&curve, MarketState { x, c: 100.0 });
        assert!(phi_hi_c <= phi_lo_c);
    }

    // Prop 3.7: positive homogeneous, Phi(tX,tC) = t*Phi(X,C).
    #[test]
    fn phi_positive_homogeneous() {
        let curve = PutPremiumCurve { k: 100.0 };
        let base = total_premium(&curve, MarketState { x: 30.0, c: 80.0 });
        let scaled = total_premium(&curve, MarketState { x: 90.0, c: 240.0 });
        assert_relative_eq!(scaled, 3.0 * base, epsilon = 1e-9);
    }
}
