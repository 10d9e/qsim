# Quantum Circuit Simulator in Rust

A **Rust-based quantum circuit simulator** supporting multiple common gates (single-qubit, multi-qubit), measurement operations (including partial measurements of subsets of qubits), and arbitrary user-defined gates. This project demonstrates the foundational ideas of quantum computing using matrix-vector multiplication on a state vector representing \(|\psi\rangle\).

> **Table of Contents**  
> - [Features](#features)  
> - [Prerequisites](#prerequisites)  
> - [Project Structure](#project-structure)  
> - [Installation and Building](#installation-and-building)  
> - [Usage](#usage)  
> - [Complete Source Code](#complete-source-code)  
>   - [Cargo.toml](#cargotoml)  
>   - [src/main.rs](#srcmainrs)  
> - [Examples](#examples)  
> - [Contributing](#contributing)  
> - [License](#license)

---

## Features

1. **N-Qubit State Representation**  
   - Stores the wavefunction \(|\psi\rangle\) in a length-\(2^N\) complex vector.

2. **Comprehensive Single-Qubit Gates**  
   - \(X, Y, Z\), **Hadamard** (H), **Phase** gates (S, S\(\dagger\), T, T\(\dagger\), etc.), **sqrt(X)** (SX).

3. **Multi-Qubit Gates**  
   - **CNOT**, **CZ**, **SWAP**, **Toffoli** (CCNOT), and an **Arbitrary** multi-qubit gate with a custom \(2^k \times 2^k\) matrix.

4. **Measurement**  
   - Single-qubit measurement (collapses the wavefunction to \(|0\rangle\) or \(|1\rangle\)).  
   - Multi-qubit measurement on a subset of qubits.  
   - Measurement of all qubits at once, collapsing the state to a single basis.

5. **Circuit Model**  
   - Construct a `QuantumCircuit` by adding gates in sequence.  
   - Apply the circuit to a `QuantumState`, which automatically updates all amplitudes.

6. **Extensibility**  
   - Easily add new gates by providing a suitable matrix.  
   - General multi-qubit gate application logic is included.

---

## Prerequisites

- [Rust](https://www.rust-lang.org/) (1.56+ recommended)  
- [Cargo](https://doc.rust-lang.org/cargo/) (typically installed alongside Rust)  
- (Optional) [Git](https://git-scm.com/) for version control and pushing to GitHub  

---

## Project Structure

```
quantum-circuit-simulator/
├── Cargo.toml
├── README.md  <-- (you are here)
└── src
    └── main.rs
```

- **Cargo.toml** – Project metadata and dependencies.  
- **src/main.rs** – Core logic of the simulator: `QuantumState`, gates, measurement functions, circuit building, etc.  

---

## Installation and Building

1. **Clone or download** this repository:
   ```bash
   git clone https://github.com/<your-username>/quantum-circuit-simulator.git
   ```
   or download and extract the ZIP file.

2. **Navigate** into the project directory:
   ```bash
   cd quantum-circuit-simulator
   ```

3. **Build** with Cargo:
   ```bash
   cargo build
   ```
   This fetches dependencies and compiles the project.

4. **Run** the simulator example:
   ```bash
   cargo run
   ```
   You’ll see printed output describing the initial state, the final state after gates, measurement results, etc.

---

## Usage

1. **Create a quantum state**:
   ```rust
   let mut qstate = QuantumState::new(3); // 3-qubit state
   ```
2. **Build a circuit** by adding gates:
   ```rust
   let mut circuit = QuantumCircuit::new(3);

// Example gates:
circuit.add_gate(QuantumGate::new_h(), &[0]);
circuit.add_gate(QuantumGate::new_x(), &[1]);
circuit.add_gate(QuantumGate::new_cnot(), &[0, 1]);
circuit.add_gate(QuantumGate::new_toffoli(), &[0, 1, 2]);
```
3. **Apply** the circuit to the state:
   ```rust
   circuit.apply(&mut qstate);
   ```
4. **Measure** qubits:
   ```rust
   // Measure a single qubit
   let outcome = qstate.measure_qubit(2);
   println!("Measured qubit 2: {}", outcome);

   // Or measure multiple qubits
   let outcome_2bits = qstate.measure_qubits(&[0, 1]);
   println!("Measured qubits [0,1]: {:#b}", outcome_2bits);

   // Or measure all qubits
   let final_measurement = qstate.measure_all();
   println!("Measured all qubits: {:#b}", final_measurement);
   ```
5. **Print** the state:
   ```rust
   qstate.print_state(); // Lists all amplitudes in binary basis
   ```

## Examples

- **Run the simulator example directly**:

  ```bash
  cargo run --example simulator
  ```

  Output will show the **initial state** (\(|000\rangle\) for 3 qubits), then the **state** after your specified gates (Hadamard, X, CNOT, etc.), and finally the results of **measurement** and the collapsed state.

- **Add new gates** (like a custom 2-qubit gate with a user-defined matrix):

  ```rust
  let custom_matrix = DMatrix::<Complex64>::from_iterator(
      4, 4,
      vec![
          Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), /* etc... */
          // fill in the 16 entries
      ]
  );
  let custom_gate = QuantumGate::new_arbitrary(custom_matrix);
  circuit.add_gate(custom_gate, &[0, 2]); // apply to qubits 0 and 2
  ```

- **Measure subsets of qubits**:

  ```rust
  let subset_result = qstate.measure_qubits(&[2, 0]);
  // This result is an integer whose bit 0 = qubit2's measurement, bit 1 = qubit0's measurement
  ```

---

## Contributing

1. **Fork** this repository on GitHub.  
2. **Create a branch** for your feature:  
   ```bash
   git checkout -b my-new-feature
   ```  
3. **Commit** your changes:  
   ```bash
   git commit -am "Add a new quantum gate or feature"
   ```  
4. **Push** your branch:  
   ```bash
   git push origin my-new-feature
   ```  
5. **Open a Pull Request** describing your changes, so it can be reviewed and merged.

---

## License

You are free to use and adapt this code in accordance with your chosen license. Commonly, open-source Rust projects are released under [MIT](https://opensource.org/licenses/MIT) or [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0). Feel free to include your own LICENSE file at the root of the repository.

---

> **Enjoy experimenting with quantum computing in Rust!** If you have questions or suggestions, feel free to open an issue or reach out.