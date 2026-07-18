pub mod black_scholes;
pub mod comparative_statics;
pub mod effective_maturity;
pub mod greeks;

pub use black_scholes::{
    exercise_boundary_call, exercise_boundary_put, price_call, price_put, AmpoParams,
};
