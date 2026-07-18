//! Property tests for the two invariants the design paper actually proves formally:
//! solvency (R_t = Phi(X_t,C_t) >= 0 always, Cor 3.12) and path independence
//! (R only depends on the current state, not the sequence of operations that got
//! there, Cor 3.11). These aren't things we're hoping are true, the paper has proofs
//! for both, so a failing test here means a bug in delta_phi or in how operations.rs
//! composes with market.rs, not a gap in the underlying math.

#[cfg(test)]
mod tests {
    use crate::market::{total_premium, MarketState};
    use crate::operations::{
        amortization_yield, buy_to_close, buy_to_open, sell_to_close, sell_to_open,
    };
    use crate::premium_curve::CallPremiumCurve;

    struct Xorshift64(u64);

    impl Xorshift64 {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        fn next_f64(&mut self, max: f64) -> f64 {
            let bits = self.next_u64();
            (bits as f64 / u64::MAX as f64) * max
        }

        fn next_range(&mut self, lo: usize, hi: usize) -> usize {
            lo + (self.next_u64() as usize) % (hi - lo)
        }
    }

    #[derive(Debug, Clone)]
    enum Action {
        Deposit(f64),
        Withdraw(f64),
        Buy(f64),
        Sell(f64),
        Amortize(f64),
    }

    fn random_action(rng: &mut Xorshift64) -> Action {
        match rng.next_range(0, 5) {
            0 => Action::Deposit(rng.next_f64(50.0)),
            1 => Action::Withdraw(rng.next_f64(50.0)),
            2 => Action::Buy(rng.next_f64(50.0)),
            3 => Action::Sell(rng.next_f64(50.0)),
            _ => Action::Amortize(rng.next_f64(5.0)),
        }
    }

    /// Applies an action to a running (state, reserves) pair, clamping the raw
    /// random amount into a valid range for the current state so we never violate
    /// X<=C or go negative. Utilization is capped short of 1.0 on buys, the call
    /// curve's premium genuinely diverges there and near-infinite premiums would
    /// swamp the final equality check with floating point noise.
    fn apply(
        curve: &CallPremiumCurve,
        state: MarketState,
        reserves: f64,
        action: &Action,
    ) -> (MarketState, f64) {
        const Q: f64 = 0.1;
        match *action {
            Action::Deposit(raw) => {
                let c = raw.max(0.0);
                let delta = sell_to_open(curve, state, c);
                (
                    MarketState {
                        x: state.x,
                        c: state.c + c,
                    },
                    reserves + delta,
                )
            }
            Action::Withdraw(raw) => {
                let max_withdraw = (state.c - state.x / 0.95).max(0.0);
                let c = raw.min(max_withdraw).max(0.0);
                let delta = buy_to_close(curve, state, c);
                (
                    MarketState {
                        x: state.x,
                        c: state.c - c,
                    },
                    reserves + delta,
                )
            }
            Action::Buy(raw) => {
                let headroom = (state.c * 0.95 - state.x).max(0.0);
                let x = raw.min(headroom).max(0.0);
                let delta = buy_to_open(curve, state, x);
                (
                    MarketState {
                        x: state.x + x,
                        c: state.c,
                    },
                    reserves + delta,
                )
            }
            Action::Sell(raw) => {
                let x = raw.min(state.x).max(0.0);
                let delta = sell_to_close(curve, state, x);
                (
                    MarketState {
                        x: state.x - x,
                        c: state.c,
                    },
                    reserves + delta,
                )
            }
            Action::Amortize(elapsed) => {
                if state.x == 0.0 {
                    return (state, reserves);
                }
                let decayed = state.x * (1.0 - (-Q * elapsed).exp());
                let delta = amortization_yield(curve, state, Q, elapsed);
                (
                    MarketState {
                        x: state.x - decayed,
                        c: state.c,
                    },
                    reserves + delta,
                )
            }
        }
    }

    fn run_trial(seed: u64) {
        let mut rng = Xorshift64(seed);
        let curve = CallPremiumCurve;
        let mut state = MarketState { x: 0.0, c: 0.0 };
        let mut reserves = 0.0_f64;

        let n_actions = rng.next_range(1, 40);
        for _ in 0..n_actions {
            let action = random_action(&mut rng);
            let (next_state, next_reserves) = apply(&curve, state, reserves, &action);
            state = next_state;
            reserves = next_reserves;
        }

        let phi_final = total_premium(&curve, state);

        assert!(
            (reserves - phi_final).abs() < 1e-6,
            "seed {}: reserves {} != Phi(final state) {}",
            seed,
            reserves,
            phi_final
        );
        assert!(
            reserves >= -1e-9,
            "seed {}: reserves went negative: {}",
            seed,
            reserves
        );
    }

    #[test]
    fn solvency_and_path_independence_hold_under_random_operation_sequences() {
        for seed in 1..500u64 {
            run_trial(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1);
        }
    }

    #[test]
    fn two_orderings_to_the_same_state_produce_the_same_reserves() {
        let curve = CallPremiumCurve;

        let s0 = MarketState { x: 0.0, c: 0.0 };
        let r_deposit_first = sell_to_open(&curve, s0, 100.0);
        let s1 = MarketState { x: 0.0, c: 100.0 };
        let r_buy_second = buy_to_open(&curve, s1, 30.0);
        let path_a_reserves = r_deposit_first + r_buy_second;

        let direct = total_premium(&curve, MarketState { x: 30.0, c: 100.0 })
            - total_premium(&curve, MarketState { x: 0.0, c: 0.0 });

        assert!((path_a_reserves - direct).abs() < 1e-9);
    }
}
