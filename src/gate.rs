use nalgebra::DMatrix;
use num_complex::Complex64;

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
