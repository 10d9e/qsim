use qsim::prelude::*;

/// Example main function demonstrating usage of the expanded simulator.
fn main() {
    // We'll do a 3-qubit system.
    let mut qstate = QuantumState::new(3);

    println!("Initial State:");
    qstate.print_state();

    // Build a circuit that does:
    // 1) H on qubit 0
    // 2) X on qubit 1
    // 3) CNOT(q0->q1)
    // 4) Toffoli(q0,q1->q2)
    // 5) Measure qubit 2
    // We'll show final state after partial measurement, etc.

    let mut circuit = QuantumCircuit::new(3);

    // 1) H on qubit 0
    circuit.add_gate(QuantumGate::new_h(), &[0]);

    // 2) X on qubit 1
    circuit.add_gate(QuantumGate::new_x(), &[1]);

    // 3) CNOT(q0->q1)
    circuit.add_gate(QuantumGate::new_cnot(), &[0, 1]);

    // 4) Toffoli(q0,q1->q2)
    circuit.add_gate(QuantumGate::new_toffoli(), &[0, 1, 2]);

    // Apply the circuit
    circuit.apply(&mut qstate);

    println!("\nState After Circuit (before measurement):");
    qstate.print_state();

    // Measure qubit 2
    let result_q2 = qstate.measure_qubit(2);
    println!("\nMeasurement of qubit 2 yielded: {}", result_q2);

    // Show state after collapse
    println!("\nState After Measuring Qubit 2:");
    qstate.print_state();

    // Let's measure the remaining qubits [0, 1] to see the final bit pattern
    let result_01 = qstate.measure_qubits(&[0, 1]);
    // result_01 is a 2-bit integer where bit0 is the measurement of qubit0 and bit1 is the measurement of qubit1
    println!("Measurement of qubits [0,1]: {:#b} (binary)", result_01);

    // Alternatively, we can measure_all:
    // (But note we've already collapsed some qubits, so the full state is basically determined).
    let all_result = qstate.measure_all();
    println!("Measurement of all qubits: {:#b}", all_result);
}
