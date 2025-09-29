#include "test_framework.hpp"
#include <OpSum.h>
#include <sstream>

using namespace pmrqmc;

bool test_opsum_construction() {
    OpSum opsum("test");
    TEST_ASSERT(opsum.get_name() == "test", "Constructor with name failed");
    TEST_ASSERT(opsum.get_terms().empty(), "New OpSum should be empty");
    TEST_ASSERT(opsum.get_max_qubit_index() == 0, "Empty OpSum should have 0 qubits");

    OpSum opsum2;
    TEST_ASSERT(opsum2.get_name().empty(), "Default constructor should have empty name");

    return true;
}

bool test_opsum_add_terms() {
    OpSum opsum;

    // Add X operator on qubit 0
    opsum.add(1.0, "X", 0);
    TEST_ASSERT(opsum.get_terms().size() == 1, "Should have 1 term after adding X_0");
    TEST_ASSERT(opsum.get_max_qubit_index() == 1, "Should have 1 qubit (0-based)");

    // Add Y operator on qubit 1
    opsum.add(2.0, "Y", 1);
    TEST_ASSERT(opsum.get_terms().size() == 2, "Should have 2 terms");
    TEST_ASSERT(opsum.get_max_qubit_index() == 2, "Should have 2 qubits");

    // Check coefficients
    auto terms = opsum.get_terms();
    TEST_ASSERT(terms[0].coefficient == std::complex<double>(1.0, 0.0), "First term should be 1.0");
    TEST_ASSERT(terms[1].coefficient == std::complex<double>(0.0, 2.0), "Y term should be i*2.0");

    return true;
}

bool test_opsum_complex_coefficients() {
    OpSum opsum;

    // Add with complex coefficients and initializer lists
    opsum.add(std::complex<double>(1.0, -1.0), {{"X", 0}, {"Z", 1}});

    auto terms = opsum.get_terms();
    TEST_ASSERT(terms.size() == 1, "Should have 1 term");
    TEST_ASSERT(terms[0].coefficient == std::complex<double>(1.0, -1.0), "Coefficient should be complex");
    TEST_ASSERT(terms[0].pauli_ops.size() == 2, "Should have 2 Pauli ops");

    // Check per qubit 2
    TEST_ASSERT(opsum.get_max_qubit_index() == 2, "Should have 2 qubits");

    return true;
}

bool test_opsum_spin_operators() {
    OpSum opsum;

    // Test Sx operator (should convert to X with 1/2 factor)
    opsum.add(2.0, "Sx", 0); // Total coefficient for Sx should be 1.0
    opsum.add(2.0, "Sy", 0); // Total coefficient for Sy should be i*1.0
    opsum.add(2.0, "Sz", 0); // Total coefficient for Sz should be 1.0

    auto terms = opsum.get_terms();
    TEST_ASSERT(terms.size() == 3, "Should have 3 terms for spin operators");

    // Check coefficients - each S operator gets multiplied by 0.5
    // Sx: 2.0 * 0.5 = 1.0 (real)
    // Sy: i * (2.0 * 0.5) = i*1.0 (imaginary)
    // Sz: 1.0 * (2.0 * 0.5) = 1.0 (real)
    TEST_ASSERT_NEAR(std::real(terms[0].coefficient), 1.0, 1e-10, "Sx coefficient should be 1.0");
    TEST_ASSERT_NEAR(std::imag(terms[0].coefficient), 0.0, 1e-10, "Sx coefficient should be real");

    TEST_ASSERT_NEAR(std::real(terms[1].coefficient), 0.0, 1e-10, "Sy real part should be 0");
    TEST_ASSERT_NEAR(std::imag(terms[1].coefficient), 1.0, 1e-10, "Sy imag part should be 1.0");

    TEST_ASSERT_NEAR(std::real(terms[2].coefficient), 1.0, 1e-10, "Sz coefficient should be 1.0");
    TEST_ASSERT_NEAR(std::imag(terms[2].coefficient), 0.0, 1e-10, "Sz coefficient should be real");

    // Check that all got converted to X, Y, Z respectively
    TEST_ASSERT(terms[0].pauli_ops[0].first == "X", "Sx should become X");
    TEST_ASSERT(terms[1].pauli_ops[0].first == "Y", "Sy should become Y");
    TEST_ASSERT(terms[2].pauli_ops[0].first == "Z", "Sz should become Z");

    return true;
}

bool test_opsum_clear() {
    OpSum opsum;
    opsum.add(1.0, "X", 0);
    opsum.add(2.0, "Y", 1);

    TEST_ASSERT(opsum.get_terms().size() == 2, "Should have 2 terms before clear");
    TEST_ASSERT(opsum.get_max_qubit_index() == 2, "Should have 2 qubits before clear");
    opsum.clear();
    TEST_ASSERT(opsum.get_terms().empty(), "Should be empty after clear");
    TEST_ASSERT(opsum.get_max_qubit_index() == 0, "Should have 0 qubits after clear");

    return true;
}

bool test_opsum_error_handling() {
    OpSum opsum;

    // Test invalid Pauli operator
    try {
        opsum.add(1.0, "W", 0); // Invalid operator
        TEST_ASSERT(false, "Should have thrown for invalid operator");
    } catch (const std::invalid_argument&) {
        // Expected
    }

    // Test negative qubit index
    try {
        opsum.add(1.0, "X", -1);
        TEST_ASSERT(false, "Should have thrown for negative qubit");
    } catch (const std::invalid_argument&) {
        // Expected
    }

    // Test empty term
    try {
        opsum.add(1.0, std::vector<std::pair<std::string, int>>{}); // Empty
        TEST_ASSERT(false, "Should have thrown for empty term");
    } catch (const std::invalid_argument&) {
        // Expected
    }

    return true;
}

bool test_opsum_template_add() {
    OpSum opsum;

    // Test template-style addition
    opsum.add(1.5, "X", 0, "Z", 1);  // X_0 * Z_1

    auto terms = opsum.get_terms();
    TEST_ASSERT(terms.size() == 1, "Should have 1 term");
    TEST_ASSERT(terms[0].coefficient == std::complex<double>(1.5, 0.0), "Coefficient should be 1.5");
    TEST_ASSERT(terms[0].pauli_ops.size() == 2, "Should have 2 Pauli operations");
    TEST_ASSERT(terms[0].pauli_ops[0] == std::make_pair(std::string("X"), 0), "First op should be X_0");
    TEST_ASSERT(terms[0].pauli_ops[1] == std::make_pair(std::string("Z"), 1), "Second op should be Z_1");

    return true;
}

bool test_opsum_print() {
    OpSum opsum("test_hamiltonian");
    opsum.add(2.0, "X", 0, "Z", 1);

    std::ostringstream oss;
    std::streambuf* old_buf = std::cout.rdbuf(oss.rdbuf());
    opsum.print();
    std::cout.rdbuf(old_buf);

    std::string output = oss.str();
    TEST_ASSERT(output.find("Hamiltonian with 1 terms:") != std::string::npos, "Output should contain term count");
    TEST_ASSERT(output.find("2") != std::string::npos, "Output should contain coefficient");
    TEST_ASSERT(output.find("X_0") != std::string::npos, "Output should contain X_0");
    TEST_ASSERT(output.find("Z_1") != std::string::npos, "Output should contain Z_1");

    return true;
}

int main() {
    TestSuite suite;

    suite.add_test("OpSum Construction", test_opsum_construction);
    suite.add_test("OpSum Add Terms", test_opsum_add_terms);
    suite.add_test("OpSum Complex Coefficients", test_opsum_complex_coefficients);
    suite.add_test("OpSum Spin Operators", test_opsum_spin_operators);
    suite.add_test("OpSum Clear", test_opsum_clear);
    suite.add_test("OpSum Error Handling", test_opsum_error_handling);
    suite.add_test("OpSum Template Add", test_opsum_template_add);
    suite.add_test("OpSum Print", test_opsum_print);

    return suite.run_all() ? 0 : 1;
}