use nalgebra::DVector;
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
        let _prob1 = 1.0 - prob0;

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
