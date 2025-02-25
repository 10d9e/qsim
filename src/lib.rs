use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use rand::Rng; // for measurement randomness

/// A QuantumState represents the state of N qubits as a vector of complex amplitudes.
///
/// The dimension of the state vector is 2^N. Indexing in the vector corresponds
/// to the binary representation of the qubit states, e.g., for N=3:
///   index = 0 -> |000>, index = 1 -> |001>, ..., index = 7 -> |111>.
///
/// We store the state as a DVector<Complex64> from the nalgebra library.
#[derive(Debug, Clone)]
pub struct QuantumState {
    /// The vector of complex amplitudes for the state of all qubits.
    pub state: DVector<Complex64>,
    /// Number of qubits.
    pub num_qubits: usize,
}

impl QuantumState {
    /// Creates a new QuantumState for N qubits in the |0...0> (all-zero) initial state.
    pub fn new(num_qubits: usize) -> Self {
        let dimension = 1 << num_qubits; // 2^num_qubits
        let mut state_vector =
            DVector::<Complex64>::from_element(dimension, Complex64::new(0.0, 0.0));
        // Amplitude of |000...0> is 1, everything else is 0.
        state_vector[0] = Complex64::new(1.0, 0.0);

        QuantumState {
            state: state_vector,
            num_qubits,
        }
    }

    /// Normalizes the quantum state's amplitudes so that the sum of squared magnitudes is 1.
    pub fn normalize(&mut self) {
        let norm_sqr: f64 = self.state.iter().map(|c| c.norm_sqr()).sum();
        if norm_sqr > 0.0 {
            let norm = norm_sqr.sqrt();
            for amp in self.state.iter_mut() {
                *amp /= Complex64::new(norm, 0.0);
            }
        }
    }

    /// Prints the amplitudes of the state in a nicely readable format.
    /// Helpful for debugging or demonstration.
    pub fn print_state(&self) {
        println!("QuantumState ({} qubits):", self.num_qubits);
        for (index, amplitude) in self.state.iter().enumerate() {
            // Print basis index in binary with leading zeros, for clarity.
            println!(
                "|{:0width$b}>: {}",
                index,
                amplitude,
                width = self.num_qubits
            );
        }
    }

    /// Measures a single qubit, returning 0 or 1, and collapses the state accordingly.
    ///
    /// # Arguments
    /// * `qubit` - The index of the qubit to measure.
    ///
    /// # Returns
    /// * `0` or `1` as the measurement result.
    pub fn measure_qubit(&mut self, qubit: usize) -> u8 {
        let dimension = 1 << self.num_qubits;

        // Probability of measuring qubit = 0 or 1
        let mut prob0 = 0.0;
        // Sum the squared magnitude of amplitudes for states where 'qubit' bit is 0.
        for i in 0..dimension {
            let bit = (i >> qubit) & 1;
            if bit == 0 {
                prob0 += self.state[i].norm_sqr();
            }
        }
        let prob1 = 1.0 - prob0;

        // Sample a random outcome from [0, 1).
        let mut rng = rand::rng();
        let sample: f64 = rng.random();

        let outcome: u8 = if sample < prob0 { 0 } else { 1 };

        // Collapse the state: zero out amplitudes of states that do not match the outcome.
        for i in 0..dimension {
            let bit = (i >> qubit) & 1;
            if bit != outcome as usize {
                self.state[i] = Complex64::new(0.0, 0.0);
            }
        }

        // Now normalize the post-measurement state.
        self.normalize();

        outcome
    }

    /// Measures multiple qubits at once, returning an integer whose bits correspond
    /// to the measurement results of those qubits in the order they appear in `qubits`.
    ///
    /// For example, if `qubits = [2, 0]` and the measurement result is qubit2=1 and qubit0=0,
    /// we return 0b10 = 2. (We place the measurement of qubit2 in the more significant bit
    /// of the returned result if it's listed first in the slice.)
    ///
    /// This function collapses the wavefunction to the subspace consistent with the measurement result.
    ///
    /// # Arguments
    /// * `qubits` - slice of qubit indices to measure
    ///
    /// # Returns
    /// * An integer with bits of the measurement results in the order of `qubits[0], qubits[1], ...`.
    pub fn measure_qubits(&mut self, qubits: &[usize]) -> u64 {
        // We'll compute probabilities for each of the 2^(qubits.len()) possible outcomes,
        // randomly select an outcome, then collapse the wavefunction accordingly.
        let num_target_qubits = qubits.len();
        if num_target_qubits == 0 {
            return 0;
        }

        let dimension = 1 << self.num_qubits;
        let num_outcomes = 1 << num_target_qubits;

        // Step 1: Compute the probability for each outcome (0..2^k-1).
        let mut probabilities = vec![0.0; num_outcomes];
        for i in 0..dimension {
            // For each basis state i, determine the bits of i for the measured qubits.
            let mut outcome_index = 0_usize;
            for (pos, qubit) in qubits.iter().enumerate() {
                let bit = ((i >> qubit) & 1) << pos;
                outcome_index |= bit;
            }
            probabilities[outcome_index] += self.state[i].norm_sqr();
        }

        // Step 2: Sample from these probabilities.
        let mut rng = rand::rng();
        let sample: f64 = rng.random();
        let mut cumulative = 0.0;
        let mut selected_outcome = 0_usize;
        for (idx, &prob) in probabilities.iter().enumerate() {
            cumulative += prob;
            if sample < cumulative {
                selected_outcome = idx;
                break;
            }
        }

        // Step 3: Collapse the state. Zero out amplitudes that do not match selected_outcome.
        for i in 0..dimension {
            // Build the outcome for the qubits in question from the index i
            let mut outcome_index = 0_usize;
            for (pos, qubit) in qubits.iter().enumerate() {
                let bit = ((i >> qubit) & 1) << pos;
                outcome_index |= bit;
            }
            if outcome_index != selected_outcome {
                self.state[i] = Complex64::new(0.0, 0.0);
            }
        }

        // Renormalize
        self.normalize();

        // Convert selected_outcome from a `usize` to `u64`
        // but keep in mind that the bit positions in selected_outcome match the order in `qubits`.
        selected_outcome as u64
    }

    /// Measures all qubits and returns the result as an integer from 0 to 2^N - 1,
    /// corresponding to the bits of each qubit in order (qubit0 is LSB).
    ///
    /// This function collapses the wavefunction to a single basis state.
    pub fn measure_all(&mut self) -> u64 {
        let dimension = 1 << self.num_qubits;

        // Probability of each basis state is amplitude.norm_sqr()
        let mut probabilities = vec![0.0; dimension];
        for i in 0..dimension {
            probabilities[i] = self.state[i].norm_sqr();
        }

        // Sample an outcome from these probabilities
        let mut rng = rand::rng();
        let sample: f64 = rng.random();
        let mut cumulative = 0.0;
        let mut selected_state = 0_usize;
        for (i, &p) in probabilities.iter().enumerate() {
            cumulative += p;
            if sample < cumulative {
                selected_state = i;
                break;
            }
        }

        // Collapse the state
        for i in 0..dimension {
            if i != selected_state {
                self.state[i] = Complex64::new(0.0, 0.0);
            }
        }
        self.normalize();

        selected_state as u64
    }
}

/// QuantumGate represents various quantum gates by storing their corresponding
/// matrix. Single-qubit gates have a 2x2 matrix, two-qubit gates have 4x4, etc.
/// For an arbitrary k-qubit gate, we store a 2^k x 2^k matrix in `Arbitrary`.
#[derive(Debug, Clone)]
pub enum QuantumGate {
    // --- Single-qubit gates ---
    X(DMatrix<Complex64>),
    Y(DMatrix<Complex64>),
    Z(DMatrix<Complex64>),
    H(DMatrix<Complex64>),
    S(DMatrix<Complex64>),
    Sdg(DMatrix<Complex64>), // S-dagger
    T(DMatrix<Complex64>),
    Tdg(DMatrix<Complex64>),   // T-dagger
    SX(DMatrix<Complex64>),    // sqrt(X)
    Phase(DMatrix<Complex64>), // a general phase shift, e.g. Rz(phi) on a single qubit

    // --- Two-qubit gates ---
    CNOT(DMatrix<Complex64>),
    CZ(DMatrix<Complex64>),
    SWAP(DMatrix<Complex64>),

    // --- Three-qubit gate (example) ---
    Toffoli(DMatrix<Complex64>), // CCNOT

    // --- Arbitrary multi-qubit gate ---
    /// The matrix must be 2^k x 2^k. The gate can be applied to k qubits in a circuit.
    Arbitrary(DMatrix<Complex64>),
}

/// Associated functions for building each gate's matrix.
impl QuantumGate {
    // --- Single-qubit gates ---

    /// Pauli-X gate
    pub fn new_x() -> Self {
        let mut mat = DMatrix::<Complex64>::zeros(2, 2);
        mat[(0, 1)] = Complex64::new(1.0, 0.0);
        mat[(1, 0)] = Complex64::new(1.0, 0.0);
        QuantumGate::X(mat)
    }

    /// Pauli-Y gate
    pub fn new_y() -> Self {
        let i = Complex64::new(0.0, 1.0);
        let mut mat = DMatrix::<Complex64>::zeros(2, 2);
        mat[(0, 1)] = -i;
        mat[(1, 0)] = i;
        QuantumGate::Y(mat)
    }

    /// Pauli-Z gate
    pub fn new_z() -> Self {
        let mut mat = DMatrix::<Complex64>::zeros(2, 2);
        mat[(0, 0)] = Complex64::new(1.0, 0.0);
        mat[(1, 1)] = Complex64::new(-1.0, 0.0);
        QuantumGate::Z(mat)
    }

    /// Hadamard gate
    pub fn new_h() -> Self {
        let val = 1.0 / (2.0_f64).sqrt();
        let mut mat = DMatrix::<Complex64>::zeros(2, 2);
        mat[(0, 0)] = Complex64::new(val, 0.0);
        mat[(0, 1)] = Complex64::new(val, 0.0);
        mat[(1, 0)] = Complex64::new(val, 0.0);
        mat[(1, 1)] = Complex64::new(-val, 0.0);
        QuantumGate::H(mat)
    }

    /// S gate: phase pi/2
    ///  [1,  0]
    ///  [0,  i]
    pub fn new_s() -> Self {
        let i = Complex64::new(0.0, 1.0);
        let mut mat = DMatrix::<Complex64>::identity(2, 2);
        mat[(1, 1)] = i;
        QuantumGate::S(mat)
    }

    /// S-dagger gate: phase -pi/2
    ///  [1,  0]
    ///  [0, -i]
    pub fn new_sdg() -> Self {
        let neg_i = Complex64::new(0.0, -1.0);
        let mut mat = DMatrix::<Complex64>::identity(2, 2);
        mat[(1, 1)] = neg_i;
        QuantumGate::Sdg(mat)
    }

    /// T gate: phase pi/4
    ///  [1, 0]
    ///  [0, exp(i pi/4)]
    pub fn new_t() -> Self {
        let exp_ipi4 = Complex64::new(
            std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
        );
        // That's e^{i pi/4} = cos(pi/4) + i sin(pi/4) = 1/sqrt(2) + i/sqrt(2).
        let mut mat = DMatrix::<Complex64>::identity(2, 2);
        mat[(1, 1)] = exp_ipi4;
        QuantumGate::T(mat)
    }

    /// T-dagger gate: phase -pi/4
    ///  [1, 0]
    ///  [0, exp(-i pi/4)]
    pub fn new_tdg() -> Self {
        let exp_minus_ipi4 = Complex64::new(
            std::f64::consts::FRAC_1_SQRT_2,
            -std::f64::consts::FRAC_1_SQRT_2,
        );
        let mut mat = DMatrix::<Complex64>::identity(2, 2);
        mat[(1, 1)] = exp_minus_ipi4;
        QuantumGate::Tdg(mat)
    }

    /// SX gate: sqrt(X) gate
    ///  1/2^(1/2) * [ 1 + i, 1 - i ]
    ///              [ 1 - i, 1 + i ]
    pub fn new_sx() -> Self {
        let factor = 1.0 / (2.0_f64).sqrt();
        let i = Complex64::new(0.0, 1.0);

        let one = Complex64::new(1.0, 0.0);
        let mut mat = DMatrix::<Complex64>::zeros(2, 2);

        // (1 + i)
        let one_plus_i = one + i;
        // (1 - i)
        let one_minus_i = one - i;

        mat[(0, 0)] = factor * one_plus_i;
        mat[(0, 1)] = factor * one_minus_i;
        mat[(1, 0)] = factor * one_minus_i;
        mat[(1, 1)] = factor * one_plus_i;

        QuantumGate::SX(mat)
    }

    /// Phase gate (single-qubit). This is effectively Rz(phi).
    ///   [1,           0      ]
    ///   [0,  exp(i * phi)    ]
    ///
    /// # Arguments
    /// * `phi` - The phase in radians.
    pub fn new_phase(phi: f64) -> Self {
        let phase = Complex64::from_polar(1.0, phi); // e^{i phi}
        let mut mat = DMatrix::<Complex64>::identity(2, 2);
        mat[(1, 1)] = phase;
        QuantumGate::Phase(mat)
    }

    // --- Two-qubit gates ---

    /// CNOT gate (control = qubit0, target = qubit1). 4x4 matrix
    ///  [1, 0, 0, 0]
    ///  [0, 1, 0, 0]
    ///  [0, 0, 0, 1]
    ///  [0, 0, 1, 0]
    pub fn new_cnot() -> Self {
        let mut mat = DMatrix::<Complex64>::zeros(4, 4);
        mat[(0, 0)] = Complex64::new(1.0, 0.0);
        mat[(1, 1)] = Complex64::new(1.0, 0.0);
        mat[(2, 3)] = Complex64::new(1.0, 0.0);
        mat[(3, 2)] = Complex64::new(1.0, 0.0);
        QuantumGate::CNOT(mat)
    }

    /// CZ gate (control = qubit0, target = qubit1). 4x4 matrix
    ///  [1, 0, 0,  0]
    ///  [0, 1, 0,  0]
    ///  [0, 0, 1,  0]
    ///  [0, 0, 0, -1]
    pub fn new_cz() -> Self {
        let mut mat = DMatrix::<Complex64>::identity(4, 4);
        mat[(3, 3)] = Complex64::new(-1.0, 0.0);
        QuantumGate::CZ(mat)
    }

    /// SWAP gate. 4x4 matrix
    ///  [1, 0, 0, 0]
    ///  [0, 0, 1, 0]
    ///  [0, 1, 0, 0]
    ///  [0, 0, 0, 1]
    pub fn new_swap() -> Self {
        let mut mat = DMatrix::<Complex64>::zeros(4, 4);
        mat[(0, 0)] = Complex64::new(1.0, 0.0);
        mat[(3, 3)] = Complex64::new(1.0, 0.0);
        mat[(1, 2)] = Complex64::new(1.0, 0.0);
        mat[(2, 1)] = Complex64::new(1.0, 0.0);
        QuantumGate::SWAP(mat)
    }

    // --- Three-qubit gates ---

    /// Toffoli gate (CCNOT). 8x8 matrix
    /// This gate flips the target qubit if both control qubits are 1.
    ///
    /// Standard layout: controls = qubits 0,1; target = qubit2. The matrix is:
    ///  [1, 0, 0, 0, 0, 0, 0, 0]
    ///  [0, 1, 0, 0, 0, 0, 0, 0]
    ///  [0, 0, 1, 0, 0, 0, 0, 0]
    ///  [0, 0, 0, 1, 0, 0, 0, 0]
    ///  [0, 0, 0, 0, 1, 0, 0, 0]
    ///  [0, 0, 0, 0, 0, 1, 0, 0]
    ///  [0, 0, 0, 0, 0, 0, 0, 1]
    ///  [0, 0, 0, 0, 0, 0, 1, 0]
    pub fn new_toffoli() -> Self {
        let mut mat = DMatrix::<Complex64>::identity(8, 8);
        // The only difference from identity is the last two basis states get swapped:
        // basis 6 <-> 7. (Binary 110 <-> 111).
        mat[(6, 6)] = Complex64::new(0.0, 0.0);
        mat[(7, 7)] = Complex64::new(0.0, 0.0);
        mat[(6, 7)] = Complex64::new(1.0, 0.0);
        mat[(7, 6)] = Complex64::new(1.0, 0.0);

        QuantumGate::Toffoli(mat)
    }

    // --- Arbitrary multi-qubit gate ---

    /// Create a gate from an arbitrary matrix. The user must ensure that
    /// `matrix` is 2^k x 2^k for some k.
    ///
    /// # Arguments
    /// * `matrix` - A DMatrix<Complex64> that is 2^k x 2^k in dimension.
    pub fn new_arbitrary(matrix: DMatrix<Complex64>) -> Self {
        QuantumGate::Arbitrary(matrix)
    }
}

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
