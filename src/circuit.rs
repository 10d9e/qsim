use crate::prelude::{QuantumGate, QuantumState};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// GateApplication: a gate + the qubits it acts on.
#[derive(Debug, Clone)]
pub struct GateApplication {
    pub gate: QuantumGate,
    pub qubits: Vec<usize>,
}

/// A QuantumCircuit is a list of gates to be applied in sequence.
#[derive(Debug, Clone)]
pub struct QuantumCircuit {
    pub num_qubits: usize,
    pub gates: Vec<GateApplication>,
}

impl QuantumCircuit {
    /// Create a new empty QuantumCircuit
    pub fn new(num_qubits: usize) -> Self {
        QuantumCircuit {
            num_qubits,
            gates: Vec::new(),
        }
    }

    /// Add a gate to the circuit with the specified qubit indices
    pub fn add_gate(&mut self, gate: QuantumGate, qubits: &[usize]) {
        self.gates.push(GateApplication {
            gate,
            qubits: qubits.to_vec(),
        });
    }

    /// Apply the entire circuit to the provided QuantumState in order.
    pub fn apply(&self, state: &mut QuantumState) {
        for gate_app in &self.gates {
            match &gate_app.gate {
                // Single-qubit gates
                QuantumGate::X(mat)
                | QuantumGate::Y(mat)
                | QuantumGate::Z(mat)
                | QuantumGate::H(mat)
                | QuantumGate::S(mat)
                | QuantumGate::Sdg(mat)
                | QuantumGate::T(mat)
                | QuantumGate::Tdg(mat)
                | QuantumGate::SX(mat)
                | QuantumGate::Phase(mat) => {
                    // These should act on exactly 1 qubit
                    let qubit_index = gate_app.qubits[0];
                    apply_single_qubit_gate(state, mat, qubit_index);
                }

                // Two-qubit gates
                QuantumGate::CNOT(mat) | QuantumGate::CZ(mat) | QuantumGate::SWAP(mat) => {
                    let q0 = gate_app.qubits[0];
                    let q1 = gate_app.qubits[1];
                    apply_two_qubit_gate(state, mat, q0, q1);
                }

                // Three-qubit gate: Toffoli
                QuantumGate::Toffoli(mat) => {
                    // Suppose the qubits are gate_app.qubits[0], [1], [2].
                    // The user must ensure they pass exactly 3 qubits.
                    let q0 = gate_app.qubits[0];
                    let q1 = gate_app.qubits[1];
                    let q2 = gate_app.qubits[2];
                    // We'll do a specialized apply for 3-qubit gates, or just use the multi-qubit approach:
                    apply_multi_qubit_gate(state, mat, &[q0, q1, q2]);
                }

                // Arbitrary multi-qubit gate
                QuantumGate::Arbitrary(mat) => {
                    let qubits = &gate_app.qubits;
                    apply_multi_qubit_gate(state, mat, qubits);
                }
            }
        }
    }
}

/// Apply a single-qubit gate to a QuantumState on a specific qubit index.
///
/// # Arguments
/// * `state` - The quantum state to modify.
/// * `gate_matrix` - A 2x2 matrix for the single-qubit gate.
/// * `qubit_index` - The index of the qubit (0-based).
fn apply_single_qubit_gate(
    state: &mut QuantumState,
    gate_matrix: &DMatrix<Complex64>,
    qubit_index: usize,
) {
    let num_qubits = state.num_qubits;
    let dimension = 1 << num_qubits;

    let mut new_state = DVector::<Complex64>::from_element(dimension, Complex64::new(0.0, 0.0));

    for i in 0..dimension {
        let bit = (i >> qubit_index) & 1;
        for row in 0..2 {
            let amp_factor = gate_matrix[(row, bit)];
            // Build new index j:
            let mut j = i;
            // Clear the bit at qubit_index in i
            if bit == 1 {
                j &= !(1 << qubit_index);
            }
            // Set if row==1
            if row == 1 {
                j |= 1 << qubit_index;
            }

            new_state[j] += state.state[i] * amp_factor;
        }
    }

    state.state = new_state;
}

/// Apply a two-qubit gate to a QuantumState on specific qubit indices.
/// `gate_matrix` is a 4x4 matrix. The qubits can be in any order, so this
/// function will handle that by extracting bits in the correct positions.
fn apply_two_qubit_gate(
    state: &mut QuantumState,
    gate_matrix: &DMatrix<Complex64>,
    qubit0: usize,
    qubit1: usize,
) {
    let num_qubits = state.num_qubits;
    let dimension = 1 << num_qubits;
    let mut new_state = DVector::<Complex64>::from_element(dimension, Complex64::new(0.0, 0.0));

    // For each basis state i, we identify the bits for qubit0, qubit1,
    // combine them into in_index in {0,1,2,3}, then map through the gate matrix
    // to find out the new amplitude contributions to out_index in {0,1,2,3}.
    for i in 0..dimension {
        let bit0 = (i >> qubit0) & 1;
        let bit1 = (i >> qubit1) & 1;

        // Combine bits in a consistent order.
        // Let's place qubit0 as the least significant bit, qubit1 as the next bit:
        let in_index = bit0 + 2 * bit1;

        // For each possible output pattern (0..3):
        for out_pattern in 0..4 {
            let amp_factor = gate_matrix[(out_pattern, in_index)];
            let out_bit0 = out_pattern & 1;
            let out_bit1 = (out_pattern >> 1) & 1;

            let mut j = i;
            // Clear the bits
            if bit0 == 1 {
                j &= !(1 << qubit0);
            }
            if bit1 == 1 {
                j &= !(1 << qubit1);
            }

            // Set bits accordingly
            if out_bit0 == 1 {
                j |= 1 << qubit0;
            }
            if out_bit1 == 1 {
                j |= 1 << qubit1;
            }

            new_state[j] += state.state[i] * amp_factor;
        }
    }
    state.state = new_state;
}

/// Apply a multi-qubit gate to a QuantumState on the specified qubits.
/// The matrix must be 2^k x 2^k, where k = qubits.len().
///
/// # Arguments
/// * `state` - The quantum state to modify.
/// * `gate_matrix` - A 2^k x 2^k matrix representing the gate.
/// * `qubits` - A slice of qubit indices. The order of these qubits determines
///              how we interpret bits in the in_index and out_index.
fn apply_multi_qubit_gate(
    state: &mut QuantumState,
    gate_matrix: &DMatrix<Complex64>,
    qubits: &[usize],
) {
    let num_qubits = state.num_qubits;
    let dimension = 1 << num_qubits;
    let k = qubits.len();
    let gate_dim = 1 << k; // 2^k

    // We'll build new_state by summing amplitude contributions.
    let mut new_state = DVector::<Complex64>::from_element(dimension, Complex64::new(0.0, 0.0));

    // For each basis state i, determine the sub-basis index in_index for the qubits of interest.
    // Then for each possible output pattern out_pattern in [0..2^k - 1], find the matrix entry
    // gate_matrix[out_pattern, in_index], and add to the appropriate new basis index j in [0..2^N - 1].
    for i in 0..dimension {
        // Extract bits for the relevant qubits from i, building in_index in [0..2^k - 1].
        let mut in_index = 0_usize;
        for (pos, &qb) in qubits.iter().enumerate() {
            let bit = ((i >> qb) & 1) << pos;
            in_index |= bit;
        }
        // out_pattern enumerates the possible new sub-states for these k qubits.
        for out_pattern in 0..gate_dim {
            let amp_factor = gate_matrix[(out_pattern, in_index)];

            // Build the new index j by clearing those bits in i and setting them to out_pattern.
            let mut j = i;
            // For each qubit, clear and then set accordingly.
            for (pos, &qb) in qubits.iter().enumerate() {
                let old_bit = (i >> qb) & 1;
                if old_bit == 1 {
                    j &= !(1 << qb);
                }
                // out_pattern's bit at position pos:
                let new_bit = (out_pattern >> pos) & 1;
                if new_bit == 1 {
                    j |= 1 << qb;
                }
            }
            new_state[j] += state.state[i] * amp_factor;
        }
    }

    state.state = new_state;
}
