//! The AmPO contract itself: notional decay and physical settlement.
//! Section 2 of the design paper (arXiv:2605.19146). This is the instrument, not the
//! market it trades in, see market.rs for the trading mechanism.

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OptionType {
    Call,
    Put,
}

#[derive(Debug, Clone, Copy)]
pub struct AmpoContract {
    pub option_type: OptionType,
    pub k: f64, // strike
    pub q: f64, // amortization rate, > 0
}

/// What each side owes on physical settlement. Direction depends on option type,
/// magnitude is symmetric: call exerciser pays cash for underlying, put exerciser
/// pays underlying for cash, both scaled by the same decayed notional.
#[derive(Debug, Clone, Copy)]
pub struct Settlement {
    pub exerciser_pays_cash: f64,
    pub exerciser_pays_underlying: f64,
    pub exerciser_receives_cash: f64,
    pub exerciser_receives_underlying: f64,
}

impl AmpoContract {
    /// N_t = N0 * e^{-q*elapsed}. Global time index, not contract age, this matters
    /// for fungibility, see the paper's remark right after Definition 2.2 in
    /// 2512.06505: two units bought at different times still decay identically from
    /// "now" if q is exogenous and shared, not from each unit's own purchase time.
    pub fn notional_at(&self, n0: f64, elapsed: f64) -> f64 {
        debug_assert!(elapsed >= 0.0);
        n0 * (-self.q * elapsed).exp()
    }

    /// Physical settlement for exercising `notional` units. Notional here should
    /// already be the decayed value from notional_at, this function doesn't apply
    /// decay itself, keeping the two concerns separate.
    pub fn settle(&self, notional: f64) -> Settlement {
        debug_assert!(notional >= 0.0);
        match self.option_type {
            OptionType::Call => Settlement {
                exerciser_pays_cash: notional * self.k,
                exerciser_pays_underlying: 0.0,
                exerciser_receives_cash: 0.0,
                exerciser_receives_underlying: notional,
            },
            OptionType::Put => Settlement {
                exerciser_pays_cash: 0.0,
                exerciser_pays_underlying: notional,
                exerciser_receives_cash: notional * self.k,
                exerciser_receives_underlying: 0.0,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn call() -> AmpoContract {
        AmpoContract {
            option_type: OptionType::Call,
            k: 100.0,
            q: 0.1,
        }
    }

    fn put() -> AmpoContract {
        AmpoContract {
            option_type: OptionType::Put,
            k: 100.0,
            q: 0.1,
        }
    }

    #[test]
    fn notional_decays_exponentially() {
        let c = call();
        assert_relative_eq!(c.notional_at(1.0, 0.0), 1.0);
        assert!(c.notional_at(1.0, 10.0) < c.notional_at(1.0, 1.0));
        assert_relative_eq!(c.notional_at(1.0, 1.0), (-0.1_f64).exp(), epsilon = 1e-12);
    }

    #[test]
    fn call_settlement_delivers_underlying_for_cash() {
        let s = call().settle(5.0);
        assert_relative_eq!(s.exerciser_pays_cash, 500.0);
        assert_relative_eq!(s.exerciser_receives_underlying, 5.0);
        assert_relative_eq!(s.exerciser_pays_underlying, 0.0);
        assert_relative_eq!(s.exerciser_receives_cash, 0.0);
    }

    #[test]
    fn put_settlement_delivers_cash_for_underlying() {
        let s = put().settle(5.0);
        assert_relative_eq!(s.exerciser_pays_underlying, 5.0);
        assert_relative_eq!(s.exerciser_receives_cash, 500.0);
        assert_relative_eq!(s.exerciser_pays_cash, 0.0);
        assert_relative_eq!(s.exerciser_receives_underlying, 0.0);
    }

    #[test]
    fn notional_never_negative_for_reasonable_horizons() {
        let c = call();
        assert!(c.notional_at(1.0, 1000.0) >= 0.0);
    }
}
