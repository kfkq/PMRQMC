#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <numeric>
#include "core/DividedDifferences.hpp"

using namespace pmrqmc::core;

// Pre-calculate factorials for the test
std::vector<double> factorials;
void precompute_factorials(int n) {
    factorials.resize(n + 1);
    factorials[0] = 1.0;
    for (int i = 1; i <= n; ++i) {
        factorials[i] = factorials[i - 1] * i;
    }
}

// Helper function for comparing ExExFloat with double
void check_close(const ExExFloat& result, double expected, double tol = 1e-12) {
    double diff = std::abs(result.to_double() - expected);
    // Add a relative tolerance check for larger numbers
    double rel_diff = expected != 0 ? diff / std::abs(expected) : diff;
    assert(diff < tol || rel_diff < tol);
}

void test_construction_and_initial_state() {
    std::cout << "Testing construction and initial state... ";
    DividedDifferences dd;
    assert(dd.empty());
    assert(dd.size() == 0);
    assert(dd.scaling_factor_s() == 1);
    std::cout << "PASS\n";
}

void test_add_single_element() {
    std::cout << "Testing adding a single element... ";
    DividedDifferences dd;
    dd.add_element(0.5);
    assert(!dd.empty());
    assert(dd.size() == 1);
    // For a single point, result is 0! * exp[z0] = exp(z0)
    check_close(dd.get_result(0), std::exp(0.5));
    std::cout << "PASS\n";
}

void test_add_multiple_elements() {
    std::cout << "Testing adding multiple elements... ";
    DividedDifferences dd;
    double z0 = 0.1, z1 = 0.2, z2 = 0.5;

    // Add first element
    dd.add_element(z0);
    check_close(dd.get_result(0), factorials[0] * std::exp(z0));

    // Add second element
    dd.add_element(z1);
    assert(dd.size() == 2);
    double expected_d1 = (std::exp(z1) - std::exp(z0)) / (z1 - z0);
    check_close(dd.get_result(0), factorials[0] * std::exp(z0)); // Should be unchanged
    check_close(dd.get_result(1), factorials[1] * expected_d1);

    // Add third element
    dd.add_element(z2);
    assert(dd.size() == 3);
    double d_z1_z2 = (std::exp(z2) - std::exp(z1)) / (z2 - z1);
    double expected_d2 = (d_z1_z2 - expected_d1) / (z2 - z0);
    check_close(dd.get_result(2), factorials[2] * expected_d2);
    
    std::cout << "PASS\n";
}

void test_remove_element() {
    std::cout << "Testing removing an element... ";
    DividedDifferences dd;
    double z0 = 0.1, z1 = 0.2, z2 = 0.5;
    dd.add_element(z0);
    dd.add_element(z1);
    dd.add_element(z2);

    assert(dd.size() == 3);

    // Remove last element (z2)
    dd.remove_element();
    assert(dd.size() == 2);
    
    // Check that the state is identical to when we only had z0 and z1
    double expected_d1 = (std::exp(z1) - std::exp(z0)) / (z1 - z0);
    check_close(dd.get_result(0), factorials[0] * std::exp(z0));
    check_close(dd.get_result(1), factorials[1] * expected_d1);

    // Remove another element (z1)
    dd.remove_element();
    assert(dd.size() == 1);
    check_close(dd.get_result(0), factorials[0] * std::exp(z0));

    // Remove the last element (z0)
    dd.remove_element();
    assert(dd.empty());
    assert(dd.size() == 0);

    std::cout << "PASS\n";
}

void test_clear() {
    std::cout << "Testing clear method... ";
    DividedDifferences dd;
    dd.add_element(0.1);
    dd.add_element(0.2);
    assert(dd.size() == 2);

    dd.clear();
    assert(dd.empty());
    assert(dd.size() == 0);
    assert(dd.scaling_factor_s() == 1); // Should reset to default

    // Should be able to add elements again
    dd.add_element(0.3);
    assert(dd.size() == 1);
    check_close(dd.get_result(0), std::exp(0.3));

    std::cout << "PASS\n";
}

void test_add_remove_sequence() {
    std::cout << "Testing add/remove sequence... ";
    DividedDifferences dd;
    
    dd.add_element(1.0);
    dd.add_element(2.0);
    dd.add_element(3.0);
    assert(dd.size() == 3);

    dd.remove_element(); // remove 3.0
    assert(dd.size() == 2);
    double expected_d01 = (std::exp(2.0) - std::exp(1.0)) / (2.0 - 1.0);
    check_close(dd.get_result(1), factorials[1] * expected_d01);

    dd.add_element(4.0); // add 4.0, sequence is {1.0, 2.0, 4.0}
    assert(dd.size() == 3);
    double d12_new = (std::exp(4.0) - std::exp(2.0)) / (4.0 - 2.0);
    double expected_d012_new = (d12_new - expected_d01) / (4.0 - 1.0);
    check_close(dd.get_result(2), factorials[2] * expected_d012_new);

    dd.remove_element();
    dd.remove_element();
    assert(dd.size() == 1);
    check_close(dd.get_result(0), std::exp(1.0));

    std::cout << "PASS\n";
}

int main() {
    precompute_factorials(20); // Precompute for tests
    std::cout << "Running DividedDifferences tests...\n\n";
    test_construction_and_initial_state();
    test_add_single_element();
    test_add_multiple_elements();
    test_remove_element();
    test_clear();
    test_add_remove_sequence();
    std::cout << "\nAll tests passed! DividedDifferences is working correctly.\n";
    return 0;
}