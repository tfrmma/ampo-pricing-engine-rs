//! Prices one AmPO call and put, prints the full risk picture (Greeks, sensitivity
//! to q, effective maturity), the sort of thing a trading desk would actually want
//! printed side by side rather than having to call each function separately.
//!
//! Run with: cargo run --release --example price_and_risk -p ampo-pricing
//! (--release matters here, effective_maturity's CRR tree and bisection are slow
//! in debug builds)

use ampo_pricing::black_scholes::{
    exercise_boundary_call, exercise_boundary_put, price_call, price_put, AmpoParams,
};
use ampo_pricing::comparative_statics::{
    dboundary_dq_call, dboundary_dq_put, dprice_dq_call, dprice_dq_put,
};
use ampo_pricing::effective_maturity::{effective_maturity_call, effective_maturity_put};
use ampo_pricing::greeks::{delta_call, delta_put, gamma_call, gamma_put, vega_call, vega_put};

fn main() {
    let p = AmpoParams {
        s0: 100.0,
        k: 100.0,
        r: 0.05,
        sigma: 0.5,
        q: 0.15,
    };
    println!(
        "AmPO, S0={} K={} r={} sigma={} q={}\n",
        p.s0, p.k, p.r, p.sigma, p.q
    );

    println!("-- fair value --");
    println!("call: {:.4}   put: {:.4}", price_call(&p), price_put(&p));

    println!("\n-- exercise boundary --");
    println!(
        "call: {:.4}   put: {:.4}",
        exercise_boundary_call(&p),
        exercise_boundary_put(&p)
    );

    println!("\n-- Greeks --");
    println!("           call        put");
    println!("delta   {:>8.4}   {:>8.4}", delta_call(&p), delta_put(&p));
    println!("gamma   {:>8.6}   {:>8.6}", gamma_call(&p), gamma_put(&p));
    println!("vega    {:>8.4}   {:>8.4}", vega_call(&p), vega_put(&p));

    println!("\n-- sensitivity to the amortization rate q itself --");
    println!("(this is a market-design parameter, not a hedging Greek)");
    println!(
        "dPrice/dq     call: {:>10.4}   put: {:>10.4}",
        dprice_dq_call(&p),
        dprice_dq_put(&p)
    );
    println!(
        "dBoundary/dq  call: {:>10.4}   put: {:>10.4}",
        dboundary_dq_call(&p),
        dboundary_dq_put(&p)
    );

    println!("\n-- effective maturity --");
    println!("(the T for a plain dated American option, no amortization, whose value happens");
    println!(" to match this AmPO's, see effective_maturity.rs docs for what this does and");
    println!(" doesn't mean, it's a value-matching reference point, not a real expiry)");
    let t_eff_call = effective_maturity_call(&p, 50.0, 500, 1e-3);
    let t_eff_put = effective_maturity_put(&p, 50.0, 500, 1e-3);
    println!("call: {:.3}y   put: {:.3}y", t_eff_call, t_eff_put);

    println!("\n-- q sweep: how fair value and boundary move across a few canonical rates --");
    println!("(Section 4.1 of the design paper suggests underwriters segment by a handful");
    println!(" of fixed q values rather than continuous maturities, similar to how a desk");
    println!(" might quote 30/60/90-day paper. This is what that segmentation looks like.)");
    println!(
        "{:>10}  {:>10}  {:>10}  {:>10}  {:>10}",
        "q", "call", "put", "call_bnd", "put_bnd"
    );
    for q in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0] {
        let pq = AmpoParams { q, ..p };
        println!(
            "{:>10.3}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}",
            q,
            price_call(&pq),
            price_put(&pq),
            exercise_boundary_call(&pq),
            exercise_boundary_put(&pq)
        );
    }
}
