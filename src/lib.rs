mod circuit;
mod gate;
mod state;

// export the public API wrapped in prelude module
pub mod prelude {
    pub use crate::circuit::QuantumCircuit;
    pub use crate::gate::QuantumGate;
    pub use crate::state::QuantumState;
}
