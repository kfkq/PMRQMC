#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>
#include <fstream>
#include <filesystem>

#include "core/Pauli.hpp"
#include "core/OperatorTerm.hpp"
#include "core/Hamiltonian.hpp"

using namespace pmrqmc::core;

void test_pauli_string_basics() {
    std::cout << "Testing PauliString basics... ";
    
    // Test default construction
    PauliString p1;
    assert(p1.is_identity());
    assert(p1.size() == 0);
    assert(p1.to_string() == "I");
    
    // Test string construction
    PauliString p2("X0 Z1 Y2");
    assert(p2.size() == 3);
    assert(p2.get_operator(0) == PauliOp::X);
    assert(p2.get_operator(1) == PauliOp::Z);
    assert(p2.get_operator(2) == PauliOp::Y);
    assert(p2.get_operator(3) == PauliOp::I);  // Not present
    
    std::string result = p2.to_string();
    // Result could be "X0 Z1 Y2" or different order
    assert(result.find("X0") != std::string::npos);
    assert(result.find("Z1") != std::string::npos);
    assert(result.find("Y2") != std::string::npos);
    
    std::cout << "PASS\n";
}

void test_pauli_string_arithmetic() {
    std::cout << "Testing PauliString arithmetic... ";
    
    PauliString p1("X0");
    PauliString p2("X0");
    PauliString p3 = p1 * p2;
    assert(p3.is_identity());  // X * X = I
    
    PauliString p4("X0");
    PauliString p5("Z1");
    PauliString p6 = p4 * p5;
    assert(p6.size() == 2);
    assert(p6.get_operator(0) == PauliOp::X);
    assert(p6.get_operator(1) == PauliOp::Z);
    
    // Test anti-commutation
    PauliString p7("Z0");  // Z on same qubit as p4
    assert(!p4.commutes_with(p7));  // X and Z on same qubit anti-commute
    assert(p4.commutes_with(p1));   // X and X commute
    
    // Test commutation on different qubits
    assert(p4.commutes_with(p5));     // X0 and Z1 commute (different qubits)
    
    std::cout << "PASS\n";
}

void test_pauli_string_application() {
    std::cout << "Testing PauliString application... ";
    
    // Test X operator (bit flip)
    PauliString x_op("X0");
    auto [state1, phase1] = x_op.apply("000");
    assert(state1 == "100");
    assert(std::abs(phase1 - 1.0) < 1e-15);
    
    // Test Z operator (phase flip)
    PauliString z_op("Z0");
    auto [state2, phase2] = z_op.apply("100");
    assert(state2 == "100");
    assert(std::abs(phase2 - (-1.0)) < 1e-15);
    
    // Test Y operator (bit flip + phase)
    PauliString y_op("Y0");
    auto [state3a, phase3a] = y_op.apply("000");
    assert(state3a == "100");
    assert(std::abs(phase3a - std::complex<double>(0, 1.0)) < 1e-15);
    
    auto [state3b, phase3b] = y_op.apply("100");
    assert(state3b == "000");
    assert(std::abs(phase3b - std::complex<double>(0, -1.0)) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_operator_term_basics() {
    std::cout << "Testing OperatorTerm basics... ";
    
    OperatorTerm term1;
    assert(term1.is_identity());
    assert(term1.is_diagonal());
    assert(term1.num_qubits() == 0);
    
    OperatorTerm term2(2.5, PauliString("X0 Z1"));
    assert(!term2.is_identity());
    assert(!term2.is_diagonal());
    assert(term2.num_qubits() == 2);
    assert(term2.acts_on_qubit(0));
    assert(term2.acts_on_qubit(1));
    assert(!term2.acts_on_qubit(2));
    assert(term2.max_qubit() == 1);
    
    OperatorTerm term3(std::complex<double>(1.0, 2.0), PauliString("Z0"));
    assert(!term3.is_identity());
    assert(term3.is_diagonal());
    
    std::cout << "PASS\n";
}

void test_operator_term_arithmetic() {
    std::cout << "Testing OperatorTerm arithmetic... ";
    
    OperatorTerm term(2.0, PauliString("X0"));
    OperatorTerm scaled = term * 3.0;
    assert(std::abs(scaled.coefficient().real() - 6.0) < 1e-15);
    
    term *= std::complex<double>(0.0, 1.0);
    assert(std::abs(term.coefficient().real()) < 1e-15);
    assert(std::abs(term.coefficient().imag() - 2.0) < 1e-15);
    
    // Test scalar multiplication from left
    OperatorTerm term2(1.0, PauliString("Z0"));
    OperatorTerm term3 = 4.0 * term2;
    assert(std::abs(term3.coefficient().real() - 4.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_hamiltonian_basics() {
    std::cout << "Testing Hamiltonian basics... ";
    
    std::vector<OperatorTerm> terms = {
        OperatorTerm(1.0, PauliString("Z0")),
        OperatorTerm(2.0, PauliString("X0")),
        OperatorTerm(-1.5, PauliString("Z1 X0"))
    };
    
    Hamiltonian h(terms);
    
    assert(h.num_qubits() == 2);  // Max qubit is 1, so 2 qubits
    assert(h.num_terms() == 3);
    assert(h.num_diagonal_terms() == 1);  // Only Z0
    assert(h.num_off_diagonal_terms() == 2);  // X0 and Z1 X0
    
    // Test adding terms
    h.add_term(OperatorTerm(0.5, PauliString("Y0")));
    assert(h.num_terms() == 4);
    
    std::cout << "PASS\n";
}

void test_hamiltonian_file_io() {
    std::cout << "Testing Hamiltonian file I/O... ";
    
    // Create a test file
    std::string test_filename = "test_hamiltonian.txt";
    std::ofstream test_file(test_filename);
    
    test_file << "# Test Hamiltonian\n";
    test_file << "1.0 Z0\n";
    test_file << "2.0 X0\n";
    test_file << "-1.5 Z1 X0\n";
    test_file << "0.5+1.0i Y1\n";
    test_file.close();
    
    // Load from file
    Hamiltonian h = Hamiltonian::from_file(test_filename);
    
    assert(h.num_qubits() == 2);
    assert(h.num_terms() == 4);
    assert(h.num_diagonal_terms() == 2);  // Z0 and Z1 are both diagonal
    assert(h.num_off_diagonal_terms() == 2);  // X0 and Y1 are off-diagonal
    
    // Test complex coefficient
    const auto& terms = h.get_terms();
    bool found_complex = false;
    for (const auto& term : terms) {
        if (std::abs(term.coefficient().imag() - 1.0) < 1e-15) {
            found_complex = true;
            break;
        }
    }
    assert(found_complex);
    
    // Cleanup
    std::filesystem::remove(test_filename);
    
    std::cout << "PASS\n";
}

void test_hamiltonian_properties() {
    std::cout << "Testing Hamiltonian properties... ";
    
    // Create a Hermitian Hamiltonian
    std::vector<OperatorTerm> terms = {
        OperatorTerm(1.0, PauliString("Z0")),      // Real term
        OperatorTerm(0.5, PauliString("X0")),      // Real term
        OperatorTerm(std::complex<double>(0.0, 1.0), PauliString("Y0")),  // Imaginary term
        OperatorTerm(std::complex<double>(0.0, -1.0), PauliString("Y0"))  // Conjugate
    };
    
    Hamiltonian h(terms);
    assert(h.is_hermitian());
    assert(!h.is_real());  // Has complex coefficients
    assert(!h.is_diagonal());  // Has off-diagonal terms
    
    assert(h.max_coefficient_magnitude() == 1.0);
    assert(h.min_coefficient_magnitude() == 0.5);
    
    // Test cleanup
    h.add_term(OperatorTerm(1e-16, PauliString("X1")));  // Very small term
    h.cleanup(1e-15);
    assert(h.num_terms() == 4);  // Small term should be removed
    
    std::cout << "PASS\n";
}

void test_hamiltonian_compatibility() {
    std::cout << "Testing legacy format compatibility... ";
    
    // Create a test file in legacy format
    std::string test_filename = "test_legacy_hamiltonian.txt";
    std::ofstream test_file(test_filename);
    
    test_file << "1.0 0 Z 1 Z\n";    // 1.0 * Z0 * Z1
    test_file << "2.0i 0 Y\n";       // 2.0i * Y0
    test_file << "0.5 0 X 1 X\n";    // 0.5 * X0 * X1
    test_file.close();
    
    // Load from legacy format
    Hamiltonian h = Hamiltonian::from_legacy_format(test_filename);
    
    assert(h.num_qubits() == 2);
    assert(h.num_terms() == 3);
    
    // Check that we can convert back
    std::string legacy_output = h.to_legacy_format();
    assert(legacy_output.find("1") != std::string::npos);
    assert(legacy_output.find("0") != std::string::npos);
    assert(legacy_output.find("Z") != std::string::npos);
    
    // Cleanup
    std::filesystem::remove(test_filename);
    
    std::cout << "PASS\n";
}

void test_performance_and_extreme_values() {
    std::cout << "Testing performance and extreme values... ";
    
    // Create a Hamiltonian with many terms
    std::vector<OperatorTerm> terms;
    for (int i = 0; i < 100; ++i) {
        double coeff = 1.0 / (i + 1);  // Decreasing coefficients
        PauliString pauli(std::string("X") + std::to_string(i % 10));  // Repeat qubits
        terms.emplace_back(coeff, pauli);
    }
    
    Hamiltonian h(terms);
    assert(h.num_terms() == 100);
    
    // Test sorting
    h.sort_terms();
    const auto& sorted_terms = h.get_terms();
    assert(sorted_terms[0].coefficient().real() > sorted_terms[1].coefficient().real());
    
    // Test very large coefficients
    OperatorTerm large_coeff_term(1e300, PauliString("Z0"));
    Hamiltonian h2({large_coeff_term});
    assert(h2.max_coefficient_magnitude() == 1e300);
    
    // Test very small coefficients
    OperatorTerm small_coeff_term(1e-100, PauliString("X0"));
    Hamiltonian h3({small_coeff_term});
    assert(h3.min_coefficient_magnitude() == 1e-100);
    
    std::cout << "PASS\n";
}

void test_pauli_string_equality_and_hashing() {
    std::cout << "Testing PauliString equality and hashing... ";
    
    PauliString p1("X0 Z1");
    PauliString p2("Z1 X0");
    PauliString p3("X0 Z1");
    
    // Test equality (order shouldn't matter for same operators)
    assert(p1 == p3);
    assert(p1 == p2);  // Different string representation but same operators
    
    // Test hashing
    assert(p1.hash() == p3.hash());  // Same operators should hash the same
    
    // Test different operators hash differently
    PauliString p4("Y0");
    PauliString p5("X0");
    assert(p4.hash() != p5.hash());
    
    std::cout << "PASS\n";
}

int main() {
    std::cout << "Running Physics Representation tests...\n\n";
    
    test_pauli_string_basics();
    test_pauli_string_arithmetic();
    test_pauli_string_application();
    test_operator_term_basics();
    test_operator_term_arithmetic();
    test_hamiltonian_basics();
    test_hamiltonian_file_io();
    test_hamiltonian_properties();
    test_hamiltonian_compatibility();
    test_performance_and_extreme_values();
    test_pauli_string_equality_and_hashing();
    
    std::cout << "\nAll Physics Representation tests passed!\n";
    std::cout << "\nPhysics representation system is working correctly.\n";
    return 0;
}