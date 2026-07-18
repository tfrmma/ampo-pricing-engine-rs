//! The six market operations, all expressed as delta_phi: change in the market's
//! total premium function when (X,C) moves by (dx,dc). Sign convention: positive
//! delta_phi means the market's reserves increase (trader pays the market), negative
//! means reserves decrease (market pays the trader a rebate). The paper phrases buy
//! operations as "pay delta_phi" and sell/rebate operations as "receive -delta_phi",
//! we just return the signed delta_phi everywhere and let the caller read the sign,
//! fewer places for a sign error to hide than juggling two conventions.

use crate::market::{total_premium, MarketState};
use crate::premium_curve::PremiumFunction;

/// Eq 3.1. The building block every operation below is defined in terms of.
pub fn delta_phi(curve: &impl PremiumFunction, before: MarketState, dx: f64, dc: f64) -> f64 {
    let after = MarketState {
        x: before.x + dx,
        c: before.c + dc,
    };
    total_premium(curve, after) - total_premium(curve, before)
}

/// Trader opens a long position of `x` notional-units. x in [0, C-X].
pub fn buy_to_open(curve: &impl PremiumFunction, state: MarketState, x: f64) -> f64 {
    debug_assert!((0.0..=state.c - state.x).contains(&x));
    delta_phi(curve, state, x, 0.0)
}

/// Trader closes `x` notional-units of a long position. x in [0, X].
pub fn sell_to_close(curve: &impl PremiumFunction, state: MarketState, x: f64) -> f64 {
    debug_assert!((0.0..=state.x).contains(&x));
    delta_phi(curve, state, -x, 0.0)
}

/// Underwriter posts `c` new collateral. First underwriter into an empty market
/// (X=0) gets zero rebate here, they're compensated later via amortization yield as
/// open interest accrues, not on entry.
pub fn sell_to_open(curve: &impl PremiumFunction, state: MarketState, c: f64) -> f64 {
    debug_assert!(c >= 0.0);
    delta_phi(curve, state, 0.0, c)
}

/// Underwriter withdraws `c` collateral. Capped at C-X, an underwriter can't pull
/// collateral that's backing live open interest.
pub fn buy_to_close(curve: &impl PremiumFunction, state: MarketState, c: f64) -> f64 {
    debug_assert!((0.0..=state.c - state.x).contains(&c));
    delta_phi(curve, state, 0.0, -c)
}

/// Holder exercises `x` notional-units. Both X and C drop by x, physical settlement
/// happens outside this function (see payoff.rs), this just captures the market-side
/// premium rebate to underwriters from the resulting drop in utilization.
pub fn exercise_yield(curve: &impl PremiumFunction, state: MarketState, x: f64) -> f64 {
    debug_assert!((0.0..=state.x).contains(&x));
    delta_phi(curve, state, -x, -x)
}

/// Passive yield accrued to underwriters from amortization decay between t0 and t,
/// on an open interest of X at t0. Remark 8 in the design paper: this only actually
/// gets computed lazily, on the next user interaction, but the yield itself is path
/// independent in time so it doesn't matter how long it's been sitting uncomputed.
pub fn amortization_yield(
    curve: &impl PremiumFunction,
    state: MarketState,
    q: f64,
    elapsed: f64,
) -> f64 {
    debug_assert!(q > 0.0 && elapsed >= 0.0);
    let decayed_fraction = 1.0 - (-q * elapsed).exp();
    delta_phi(curve, state, -decayed_fraction * state.x, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::premium_curve::CallPremiumCurve;

    #[test]
    fn buy_to_open_costs_the_trader() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 10.0, c: 100.0 };
        assert!(buy_to_open(&curve, state, 20.0) > 0.0);
    }

    #[test]
    fn sell_to_close_rebates_the_trader() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 30.0, c: 100.0 };
        assert!(sell_to_close(&curve, state, 20.0) < 0.0);
    }

    #[test]
    fn sell_to_open_gives_zero_rebate_on_empty_market() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 0.0, c: 0.0 };
        assert_eq!(sell_to_open(&curve, state, 50.0), 0.0);
    }

    #[test]
    fn sell_to_open_gives_positive_rebate_when_open_interest_exists() {
        // adding collateral while X>0 lowers utilization, which lowers phi, so
        // delta_phi is negative here, meaning the underwriter is compensated.
        let curve = CallPremiumCurve;
        let state = MarketState { x: 30.0, c: 100.0 };
        assert!(sell_to_open(&curve, state, 20.0) < 0.0);
    }

    #[test]
    fn round_trip_buy_then_sell_nets_to_zero() {
        // Remark 9: no arbitrage from round-tripping. Buy x then immediately sell x
        // back nets to exactly zero, no fees modeled here so this should be exact,
        // not just approximately zero.
        let curve = CallPremiumCurve;
        let state = MarketState { x: 10.0, c: 100.0 };
        let cost = buy_to_open(&curve, state, 15.0);
        let after_buy = MarketState {
            x: state.x + 15.0,
            c: state.c,
        };
        let rebate = sell_to_close(&curve, after_buy, 15.0);
        assert!((cost + rebate).abs() < 1e-9);
    }

    #[test]
    fn amortization_yield_is_zero_at_zero_elapsed_time() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 30.0, c: 100.0 };
        assert_eq!(amortization_yield(&curve, state, 0.1, 0.0), 0.0);
    }

    #[test]
    fn amortization_yield_grows_with_elapsed_time() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 30.0, c: 100.0 };
        let short = amortization_yield(&curve, state, 0.1, 1.0).abs();
        let long = amortization_yield(&curve, state, 0.1, 10.0).abs();
        assert!(long > short);
    }

    #[test]
    fn exercise_reduces_utilization_and_pays_a_rebate() {
        let curve = CallPremiumCurve;
        let state = MarketState { x: 50.0, c: 100.0 };
        assert!(exercise_yield(&curve, state, 20.0) < 0.0);
    }
}
