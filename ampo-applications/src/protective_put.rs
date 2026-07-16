//! Protective put collateral, Section 4.1. A borrower posts E units of an asset plus
//! alpha*E put options (strike K, amortization q, alpha > 1) instead of relying on
//! an exogenous, oracle-dependent LTV haircut. The put's decaying notional becomes
//! the collateral valuation, so liquidation timing is fully endogenous and
//! oracle-free.
//!
//! Note on scope: the paper gives exact formulas for collateral value and protection
//! decay, and describes Figure 1 (green/yellow/red liquidation regions) narratively,
//! but doesn't give a formula for where "safe liquidation" ends and "bad debt"
//! begins, that split depends on liquidator speed and slippage assumptions the paper
//! doesn't specify. This module implements the two things that ARE specified:
//! collateral value over time, and whether a given state is liquidatable. It does
//! NOT attempt to reconstruct the yellow/red split, that would be inventing a model
//! the paper doesn't give us.

/// A protective put collateral position: E units of the risky asset, alpha*E units
/// of notional in put options struck at K, amortizing at rate q.
#[derive(Debug, Clone, Copy)]
pub struct ProtectivePutCollateral {
    pub e: f64,     // units of the risky asset (e.g. WETH) posted
    pub alpha: f64, // over-collateralization multiplier on option notional, > 1
    pub k: f64,     // strike, USDC per unit of underlying
    pub q: f64,     // amortization rate
}

impl ProtectivePutCollateral {
    /// t -> alpha*e^{-qt}*E, the raw decaying put notional (Section 4.1, uncapped).
    pub fn protection_notional(&self, t: f64) -> f64 {
        debug_assert!(t >= 0.0);
        self.alpha * (-self.q * t).exp() * self.e
    }

    /// Collateral value in USDC: E*K*min(alpha*e^{-qt}, 1). The min() reflects that
    /// you can never exercise puts on more of the underlying than you actually hold
    /// (E units), the option protects at most 1:1 against the posted asset even
    /// while alpha*e^{-qt} > 1 early in the position's life.
    pub fn collateral_value(&self, t: f64) -> f64 {
        debug_assert!(t >= 0.0);
        self.e * self.k * (self.alpha * (-self.q * t).exp()).min(1.0)
    }

    /// The time at which the option notional decays through 1.0, i.e. when the min()
    /// in collateral_value starts binding and full 1:1 protection is no longer
    /// guaranteed. Before this, the position is fully hedged; after, the collateral
    /// value decays exponentially with the option.
    pub fn full_protection_expiry(&self) -> f64 {
        debug_assert!(self.alpha > 1.0, "alpha <= 1 means the position is never over-protected");
        self.alpha.ln() / self.q
    }
}

/// Loan-to-value at time t, given a debt balance that's grown from debt0 at
/// continuously compounded rate debt_growth_rate. LTV = debt_t / collateral_value_t.
pub fn loan_to_value(position: &ProtectivePutCollateral, debt0: f64, debt_growth_rate: f64, t: f64) -> f64 {
    let debt_t = debt0 * (debt_growth_rate * t).exp();
    let collateral_t = position.collateral_value(t);
    debt_t / collateral_t
}

/// Whether the position should be liquidated at time t, i.e. LTV has crossed the
/// governance-set threshold. This is the binary green/red check the paper's math
/// actually supports, no yellow band, see module docs.
pub fn is_liquidatable(
    position: &ProtectivePutCollateral,
    debt0: f64,
    debt_growth_rate: f64,
    liquidation_ltv: f64,
    t: f64,
) -> bool {
    loan_to_value(position, debt0, debt_growth_rate, t) >= liquidation_ltv
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn position() -> ProtectivePutCollateral {
        ProtectivePutCollateral { e: 1.0, alpha: 1.2, k: 3000.0, q: 0.003 } // 30bps/day-ish
    }

    #[test]
    fn collateral_value_starts_at_full_protection() {
        let p = position();
        // at t=0, alpha=1.2 > 1, min() caps it at 1.0, so full E*K value.
        assert_relative_eq!(p.collateral_value(0.0), p.e * p.k, epsilon = 1e-9);
    }

    #[test]
    fn collateral_value_decays_after_full_protection_expiry() {
        let p = position();
        let t_expiry = p.full_protection_expiry();
        assert!(t_expiry > 0.0);
        // just before expiry, still fully protected.
        assert_relative_eq!(p.collateral_value(t_expiry - 0.5), p.e * p.k, epsilon = 1e-6);
        // well after, decaying below full value.
        assert!(p.collateral_value(t_expiry + 50.0) < p.e * p.k);
    }

    #[test]
    fn collateral_value_never_negative() {
        let p = position();
        assert!(p.collateral_value(10000.0) >= 0.0);
    }

    #[test]
    fn loan_to_value_increases_as_debt_compounds_and_collateral_decays() {
        let p = position();
        let debt0 = 2000.0;
        let rate = 0.05;
        let ltv_early = loan_to_value(&p, debt0, rate, 10.0);
        let ltv_late = loan_to_value(&p, debt0, rate, 500.0);
        assert!(ltv_late > ltv_early);
    }

    #[test]
    fn liquidation_triggers_once_ltv_crosses_threshold() {
        let p = position();
        let debt0 = 2000.0;
        let rate = 0.05;
        let threshold = 0.83; // matches the Aave v3 WETH figure cited in the paper
        assert!(!is_liquidatable(&p, debt0, rate, threshold, 0.0));
        assert!(is_liquidatable(&p, debt0, rate, threshold, 2000.0));
    }
}
