//! De-peg insurance, Section 4.2. An at-the-peg put AmPO market (underlying: the
//! stablecoin, strike K=1, numeraire: the reserve asset) turns a discretionary PSM
//! (peg stability module) into an explicit, priced option. Holding and exercising
//! the put is exactly the DAI-for-USDC conversion a PSM offers, except now
//! underwriters get paid via amortization yield instead of the DAO eating the cost
//! for free, and third parties can underwrite it too instead of only the DAO.
//!
//! This is a thin domain wrapper, there's no new math here beyond Section 3's
//! market mechanism applied to K=1 puts. Resist the urge to add pricing logic in
//! this file, ampo-pricing already owns that.

use ampo_core::market::MarketState;
use ampo_core::payoff::{AmpoContract, OptionType, Settlement};
use ampo_core::premium_curve::PremiumFunction;

/// A de-peg insurance market is just a put AmPO struck at parity. K need not be
/// exactly 1.0 in principle (e.g. insuring against a small permanent de-peg band),
/// but 1.0 is the canonical PSM-equivalent case from the paper.
pub fn psm_put_contract(q: f64) -> AmpoContract {
    AmpoContract {
        option_type: OptionType::Put,
        k: 1.0,
        q,
    }
}

/// Converting `stablecoin_amount` units of the pegged asset into the reserve asset
/// via exercise, exactly the PSM conversion, just routed through physical
/// settlement instead of a discretionary DAO-run swap. Caller is responsible for
/// confirming the user actually holds this much decayed put notional first, this
/// function only computes the settlement legs.
pub fn convert_via_exercise(contract: &AmpoContract, stablecoin_amount: f64) -> Settlement {
    debug_assert_eq!(contract.option_type, OptionType::Put);
    contract.settle(stablecoin_amount)
}

/// The "fear index": current marginal premium P(U) at the prevailing utilization.
/// This is the paper's on-chain, oracle-free signal of implied de-peg risk, rising
/// as underwriting collateral gets scarcer relative to insurance demand.
pub fn fear_index(curve: &impl PremiumFunction, state: MarketState) -> f64 {
    debug_assert!(
        state.c > 0.0,
        "fear index undefined with no underwriting collateral"
    );
    curve.premium(state.x / state.c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ampo_core::premium_curve::PutPremiumCurve;
    use approx::assert_relative_eq;

    #[test]
    fn psm_contract_is_a_put_struck_at_parity() {
        let c = psm_put_contract(0.05);
        assert_eq!(c.option_type, OptionType::Put);
        assert_relative_eq!(c.k, 1.0);
    }

    #[test]
    fn conversion_delivers_reserve_asset_at_parity() {
        let c = psm_put_contract(0.05);
        let s = convert_via_exercise(&c, 1_000_000.0);
        // put settlement: exerciser pays underlying (the depegging stablecoin),
        // receives cash (the reserve asset), at strike parity.
        assert_relative_eq!(s.exerciser_pays_underlying, 1_000_000.0);
        assert_relative_eq!(s.exerciser_receives_cash, 1_000_000.0); // K=1.0
    }

    #[test]
    fn fear_index_rises_with_utilization() {
        let curve = PutPremiumCurve { k: 1.0 };
        let calm = fear_index(&curve, MarketState { x: 10.0, c: 1000.0 });
        let stressed = fear_index(
            &curve,
            MarketState {
                x: 900.0,
                c: 1000.0,
            },
        );
        assert!(stressed > calm);
    }

    #[test]
    fn fear_index_bounded_by_strike_as_utilization_saturates() {
        // put premium function is capped at K per Definition 3.1, so the fear index
        // can spike but never exceeds parity value, unlike the call side which
        // diverges to infinity.
        let curve = PutPremiumCurve { k: 1.0 };
        let near_full = fear_index(
            &curve,
            MarketState {
                x: 999.0,
                c: 1000.0,
            },
        );
        assert!(near_full <= 1.0);
    }
}
