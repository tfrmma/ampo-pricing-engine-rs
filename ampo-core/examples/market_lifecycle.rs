//! Walks one AmPO call market through its full lifecycle: an underwriter posts
//! collateral, a trader opens a position, time passes and amortization yield
//! accrues, the trader exercises, the underwriter withdraws. Prints the running
//! reserves R alongside Phi(X,C) at every step, live, the same invariant
//! invariants.rs checks with a fuzzer, here it's one readable narrative instead.
//!
//! Also the one place in this workspace that deliberately crosses the ampo-core /
//! ampo-pricing boundary: the market's utilization-based premium P(U) here is a
//! completely different number from the Black-Scholes fair value in ampo-pricing,
//! and this example prints both side by side on purpose, so the distinction
//! doesn't stay implicit. ampo-core itself has no dependency on ampo-pricing (see
//! the repo README), this is a dev-dependency scoped to this example only.
//!
//! Run with: cargo run --example market_lifecycle -p ampo-core

use ampo_core::market::{total_premium, MarketState};
use ampo_core::operations::{
    amortization_yield, buy_to_close, buy_to_open, exercise_yield, sell_to_open,
};
use ampo_core::payoff::{AmpoContract, OptionType};
use ampo_core::premium_curve::{CallPremiumCurve, PremiumFunction};
use ampo_pricing::{price_call, AmpoParams};

fn main() {
    let curve = CallPremiumCurve;
    let contract = AmpoContract {
        option_type: OptionType::Call,
        k: 100.0,
        q: 0.1,
    };

    // running totals, tracked exactly like invariants.rs does, but printed instead
    // of asserted
    let mut state = MarketState { x: 0.0, c: 0.0 };
    let mut reserves = 0.0_f64;

    let pricing_params = AmpoParams {
        s0: 100.0,
        k: 100.0,
        r: 0.05,
        sigma: 0.5,
        q: 0.1,
    };
    println!(
        "Black-Scholes fair value of this call (ampo-pricing): {:.4}\n",
        price_call(&pricing_params)
    );

    println!("-- underwriter posts 1000 units of collateral --");
    let delta = sell_to_open(&curve, state, 1000.0);
    state.c += 1000.0;
    reserves += delta;
    print_state("after deposit", &curve, state, reserves);

    println!("\n-- trader opens 300 units of notional --");
    let delta = buy_to_open(&curve, state, 300.0);
    state.x += 300.0;
    reserves += delta;
    println!("trader pays: {:.4}", delta);
    print_state("after open", &curve, state, reserves);
    println!(
        "market's instantaneous marginal premium P(U={:.2}): {:.6} (this is NOT the Black-Scholes fair value above, \
         it's the cost of the next marginal unit of notional given current utilization, a different concept)",
        state.x / state.c,
        curve.premium(state.x / state.c)
    );

    println!("\n-- 45 days pass, amortization yield accrues to underwriters --");
    let elapsed_years = 45.0 / 365.0;
    let decayed_notional = state.x * (1.0 - (-contract.q * elapsed_years).exp());
    let delta = amortization_yield(&curve, state, contract.q, elapsed_years);
    state.x -= decayed_notional;
    reserves += delta;
    println!(
        "open interest decayed from {:.4} to {:.4} units",
        state.x + decayed_notional,
        state.x
    );
    print_state("after amortization", &curve, state, reserves);

    println!("\n-- trader exercises half of their remaining position --");
    let exercise_amount = state.x / 2.0;
    let settlement = contract.settle(exercise_amount);
    println!(
        "physical settlement: trader pays {:.4} USDC, receives {:.4} units of underlying",
        settlement.exerciser_pays_cash, settlement.exerciser_receives_underlying
    );
    let delta = exercise_yield(&curve, state, exercise_amount);
    state.x -= exercise_amount;
    state.c -= exercise_amount;
    reserves += delta;
    print_state("after exercise", &curve, state, reserves);

    println!("\n-- underwriter withdraws all now-unused collateral --");
    let withdrawable = state.c - state.x;
    println!(
        "note: withdrawing the full {:.4} would push utilization to exactly 1.0, where Phi \
         diverges for calls by design (Definition 3.1's limit, see \
         call_premium_is_infinite_at_full_utilization in premium_curve.rs). Leaving a 1% buffer \
         instead, same caution invariants.rs's fuzz test uses.",
        withdrawable
    );
    let withdrawal = withdrawable * 0.99;
    let delta = buy_to_close(&curve, state, withdrawal);
    state.c -= withdrawal;
    reserves += delta;
    println!("underwriter receives: {:.4}", -delta);
    print_state("final", &curve, state, reserves);

    println!("\n-- invariant check --");
    let phi_final = total_premium(&curve, state);
    println!("cumulative reserves R:        {:.8}", reserves);
    println!("Phi(X_final, C_final):        {:.8}", phi_final);
    println!(
        "match (path independence):    {}",
        (reserves - phi_final).abs() < 1e-9
    );
    println!("solvent (R >= 0):              {}", reserves >= 0.0);
}

fn print_state(label: &str, curve: &CallPremiumCurve, state: MarketState, reserves: f64) {
    println!(
        "[{}]  X={:.4}  C={:.4}  U={:.4}  Phi(X,C)={:.4}  R={:.4}",
        label,
        state.x,
        state.c,
        if state.c > 0.0 {
            state.x / state.c
        } else {
            0.0
        },
        total_premium(curve, state),
        reserves
    );
}
