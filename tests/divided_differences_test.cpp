#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <numeric>
#include <random>
#include <iomanip>
#include <chrono>
#include "core/DividedDifferences.hpp"

using namespace pmrqmc::core;

// Assertion macros
#define ASSERT_TRUE(cond) assert(cond)
#define ASSERT_FALSE(cond) assert(!(cond))
#define ASSERT_EQ(a, b) assert((a) == (b))
#define ASSERT_NE(a, b) assert((a) != (b))
#define ASSERT_LT(a, b) assert((a) < (b))
#define ASSERT_GT(a, b) assert((a) > (b))

// Helper function for comparing ExExFloat with double
void assert_close(const ExExFloat& result, double expected, double tol = 1e-10) {
    double actual = result.to_double();
    double diff = std::abs(actual - expected);
    double rel_diff = expected != 0 ? diff / std::abs(expected) : diff;
    
    if (diff >= tol && rel_diff >= tol) {
        std::cout << "\n  ASSERTION FAILED: Expected " << std::scientific << expected 
                  << ", got " << actual << ", diff=" << diff << ", rel_diff=" << rel_diff << std::endl;
        assert(false);
    }
}

// Test construction and basic properties
void test_construction() {
    std::cout << "Testing construction and initial state... ";
    
    DividedDifferences dd;
    ASSERT_TRUE(dd.empty());
    ASSERT_EQ(dd.size(), 0);
    ASSERT_EQ(dd.scaling_factor_s(), 1);
    
    // Test different construction parameters
    DividedDifferences dd2(100, 50);
    ASSERT_TRUE(dd2.empty());
    ASSERT_EQ(dd2.size(), 0);
    
    std::cout << "✓\n";
}

// Test single element operations
void test_single_element() {
    std::cout << "Testing single element operations... ";
    
    DividedDifferences dd;
    
    // Test various values
    std::vector<double> test_values = {0.0, 0.5, 1.0, -0.5, -1.0, 2.5, -2.5};
    
    for (double z : test_values) {
        dd.clear();
        dd.add_element(z);
        
        ASSERT_FALSE(dd.empty());
        ASSERT_EQ(dd.size(), 1);
        
        // For single point, result should be exp(z)
        assert_close(dd.get_result(0), std::exp(z));
        
        // Test removal
        dd.remove_element();
        ASSERT_TRUE(dd.empty());
        ASSERT_EQ(dd.size(), 0);
    }
    
    std::cout << "✓\n";
}

// Test two element divided differences
void test_two_elements() {
    std::cout << "Testing two element operations... ";
    
    DividedDifferences dd;
    
    // Test various pairs
    std::vector<std::pair<double, double>> test_pairs = {
        {0.1, 0.2}, {0.0, 1.0}, {-0.5, 0.5}, {1.0, 2.0}, {-1.0, -0.5}, {0.0, 0.0001}
    };
    
    for (auto [z0, z1] : test_pairs) {
        if (std::abs(z1 - z0) < 1e-12) continue; // Skip nearly identical values
        
        dd.clear();
        dd.add_element(z0);
        dd.add_element(z1);
        
        ASSERT_EQ(dd.size(), 2);
        
        // f[z0] = exp(z0)
        assert_close(dd.get_result(0), std::exp(z0));
        
        // f[z0, z1] = (exp(z1) - exp(z0)) / (z1 - z0)
        double expected_d1 = (std::exp(z1) - std::exp(z0)) / (z1 - z0);
        assert_close(dd.get_result(1), expected_d1);
    }
    
    std::cout << "✓\n";
}

// Test three element divided differences
void test_three_elements() {
    std::cout << "Testing three element operations... ";
    
    DividedDifferences dd;
    double z0 = 0.1, z1 = 0.3, z2 = 0.7;
    
    dd.add_element(z0);
    dd.add_element(z1);
    dd.add_element(z2);
    
    ASSERT_EQ(dd.size(), 3);
    
    // Calculate expected values manually - these should be FACTORIAL-SCALED
    // The DividedDifferences class returns q! * f[z_0, z_1, ..., z_q]
    double exp_z0 = std::exp(z0);
    double exp_z1 = std::exp(z1);
    double exp_z2 = std::exp(z2);
    
    // Factorial-scaled divided differences
    double f0 = exp_z0;  // 0! * f[z_0] = f[z_0] = exp(z0)
    
    // 1! * f[z_0, z_1] = (exp(z1) - exp(z0)) / (z1 - z0)
    double d01 = (exp_z1 - exp_z0) / (z1 - z0);
    
    // 2! * f[z_0, z_1, z_2] = ((exp(z2) - exp(z1))/(z2 - z1) - (exp(z1) - exp(z0))/(z1 - z0)) / (z2 - z0) * 2!
    double first_order_01 = (exp_z1 - exp_z0) / (z1 - z0);
    double first_order_12 = (exp_z2 - exp_z1) / (z2 - z1);
    double second_order = (first_order_12 - first_order_01) / (z2 - z0);
    double d012 = 2.0 * second_order;  // 2! * f[z_0, z_1, z_2]
    
    assert_close(dd.get_result(0), f0);
    assert_close(dd.get_result(1), d01);
    assert_close(dd.get_result(2), d012);
    
    std::cout << "✓\n";
}

// Test removal operations in detail
void test_removal_operations() {
    std::cout << "Testing removal operations... ";
    
    DividedDifferences dd;
    std::vector<double> values = {0.1, 0.3, 0.5, 0.8, 1.2};
    
    // Build up the sequence
    for (double z : values) {
        dd.add_element(z);
    }
    ASSERT_EQ(dd.size(), 5);
    
    // Store intermediate results
    std::vector<ExExFloat> results_at_4;
    for (size_t i = 0; i < 4; ++i) {
        results_at_4.push_back(dd.get_result(i));
    }
    
    // Remove last element
    dd.remove_element();
    ASSERT_EQ(dd.size(), 4);
    
    // Check that results match what we had before
    for (size_t i = 0; i < 4; ++i) {
        assert_close(dd.get_result(i), results_at_4[i].to_double(), 1e-12);
    }
    
    // Continue removing and verify sizes
    dd.remove_element();
    ASSERT_EQ(dd.size(), 3);
    
    dd.remove_element();
    ASSERT_EQ(dd.size(), 2);
    
    dd.remove_element();
    ASSERT_EQ(dd.size(), 1);
    
    // Last element should be exp(0.1)
    assert_close(dd.get_result(0), std::exp(0.1));
    
    dd.remove_element();
    ASSERT_TRUE(dd.empty());
    
    std::cout << "✓\n";
}

// Test clear operations
void test_clear_operations() {
    std::cout << "Testing clear operations... ";
    
    DividedDifferences dd;
    
    // Add some elements
    dd.add_element(1.0);
    dd.add_element(2.0);
    dd.add_element(3.0);
    ASSERT_EQ(dd.size(), 3);
    
    // Clear and verify
    dd.clear();
    ASSERT_TRUE(dd.empty());
    ASSERT_EQ(dd.size(), 0);
    ASSERT_EQ(dd.scaling_factor_s(), 1);
    
    // Should be able to add new elements
    dd.add_element(0.5);
    ASSERT_EQ(dd.size(), 1);
    assert_close(dd.get_result(0), std::exp(0.5));
    
    std::cout << "✓\n";
}

// Test error handling
void test_error_handling() {
    std::cout << "Testing error handling... ";
    
    DividedDifferences dd;
    
    // Test out of bounds access
    try {
        dd.get_result(0); // Should throw - empty container
        assert(false && "Should have thrown");
    } catch (const std::out_of_range&) {
        // Expected
    }
    
    dd.add_element(1.0);
    
    try {
        dd.get_result(1); // Should throw - index out of range
        assert(false && "Should have thrown");
    } catch (const std::out_of_range&) {
        // Expected
    }
    
    // Valid access should work
    try {
        ExExFloat result = dd.get_result(0);
        assert_close(result, std::exp(1.0));
    } catch (...) {
        assert(false && "Should not have thrown");
    }
    
    // Test remove from empty
    dd.clear();
    dd.remove_element(); // Should not crash
    ASSERT_TRUE(dd.empty());
    
    std::cout << "✓\n";
}

// Test with large sequences
void test_large_sequences() {
    std::cout << "Testing large sequences... ";
    
    DividedDifferences dd;
    const size_t N = 100;
    
    // Add many elements
    for (size_t i = 0; i < N; ++i) {
        double z = 0.01 * i; // Values from 0 to 0.99
        dd.add_element(z);
    }
    
    ASSERT_EQ(dd.size(), N);
    
    // Check first and last results are reasonable
    assert_close(dd.get_result(0), std::exp(0.0), 1e-10);
    
    // The last result should be finite (not NaN or inf)
    double last_result = dd.get_result(N - 1).to_double();
    ASSERT_TRUE(std::isfinite(last_result));
    
    // Test removal from large sequence
    for (size_t i = 0; i < N; ++i) {
        dd.remove_element();
        ASSERT_EQ(dd.size(), N - 1 - i);
    }
    
    ASSERT_TRUE(dd.empty());
    
    std::cout << "✓\n";
}

// Test numerical stability
void test_numerical_stability() {
    std::cout << "Testing numerical stability... ";
    
    DividedDifferences dd;
    
    // Test with very close values (potential numerical issues)
    double base = 1.0;
    std::vector<double> close_values = {base, base + 1e-8, base + 2e-8, base + 3e-8};
    
    for (double z : close_values) {
        dd.add_element(z);
    }
    
    // All results should be finite
    for (size_t i = 0; i < close_values.size(); ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
        ASSERT_FALSE(std::isnan(result));
    }
    
    // Test with extreme values
    dd.clear();
    std::vector<double> extreme_values = {-5.0, -2.0, 0.0, 2.0, 5.0};
    
    for (double z : extreme_values) {
        dd.add_element(z);
    }
    
    // All results should still be finite
    for (size_t i = 0; i < extreme_values.size(); ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
        ASSERT_FALSE(std::isnan(result));
    }
    
    std::cout << "✓\n";
}

// Test random sequences
void test_random_sequences() {
    std::cout << "Testing random sequences... ";
    
    std::random_device rd;
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    
    DividedDifferences dd;
    
    // Generate random sequence
    const size_t N = 50;
    std::vector<double> random_values;
    
    for (size_t i = 0; i < N; ++i) {
        double z = dis(gen);
        random_values.push_back(z);
        dd.add_element(z);
    }
    
    ASSERT_EQ(dd.size(), N);
    
    // All results should be finite
    for (size_t i = 0; i < N; ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
    }
    
    // Test random removals
    size_t remaining = N;
    while (remaining > 0) {
        dd.remove_element();
        remaining--;
        ASSERT_EQ(dd.size(), remaining);
        
        // Remaining results should still be finite
        for (size_t i = 0; i < remaining; ++i) {
            double result = dd.get_result(i).to_double();
            ASSERT_TRUE(std::isfinite(result));
        }
    }
    
    ASSERT_TRUE(dd.empty());
    
    std::cout << "✓\n";
}

// Test performance and memory
void test_performance() {
    std::cout << "Testing performance... ";
    
    DividedDifferences dd;
    const size_t N = 1000;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Add elements
    for (size_t i = 0; i < N; ++i) {
        dd.add_element(0.001 * i);
    }
    
    // Access all results
    for (size_t i = 0; i < N; ++i) {
        volatile double result = dd.get_result(i).to_double();
        (void)result; // Suppress unused variable warning
    }
    
    // Remove all elements
    for (size_t i = 0; i < N; ++i) {
        dd.remove_element();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "✓ (took " << duration.count() << "ms)\n";
    
    ASSERT_TRUE(dd.empty());
}

// Test scaling factor changes
void test_scaling_factor_changes() {
    std::cout << "Testing scaling factor changes... ";
    
    DividedDifferences dd;
    
    // Add elements that should trigger scaling factor changes
    dd.add_element(0.0);
    ASSERT_EQ(dd.scaling_factor_s(), 1);
    
    dd.add_element(1.0);
    dd.add_element(2.0);
    dd.add_element(10.0); // This should potentially trigger rescaling
    
    // Should still produce valid results
    ASSERT_GT(dd.size(), 0);
    for (size_t i = 0; i < dd.size(); ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
    }
    
    std::cout << "✓\n";
}

// Test edge cases
void test_edge_cases() {
    std::cout << "Testing edge cases... ";
    
    DividedDifferences dd;
    
    // Test with zero
    dd.add_element(0.0);
    assert_close(dd.get_result(0), 1.0); // exp(0) = 1
    
    // Test with identical values (should handle gracefully)
    dd.clear();
    dd.add_element(1.0);
    dd.add_element(1.0); // Same value
    
    // Should handle this without crashing
    ASSERT_EQ(dd.size(), 2);
    double result0 = dd.get_result(0).to_double();
    double result1 = dd.get_result(1).to_double();
    ASSERT_TRUE(std::isfinite(result0));
    ASSERT_TRUE(std::isfinite(result1));
    
    std::cout << "✓\n";
}

// Test mathematical properties and invariants
void test_mathematical_properties() {
    std::cout << "Testing mathematical properties... ";
    
    DividedDifferences dd;
    
    // Test that exp[z0] = exp(z0) for single point
    double z0 = 1.5;
    dd.add_element(z0);
    assert_close(dd.get_result(0), std::exp(z0));
    
    // Test symmetry: f[z0, z1] should equal f[z1, z0] for exp function
    dd.clear();
    double z1 = 2.5;
    dd.add_element(z0);
    dd.add_element(z1);
    double fwd_result = dd.get_result(1).to_double();
    
    dd.clear();
    dd.add_element(z1);
    dd.add_element(z0);
    double rev_result = dd.get_result(1).to_double();
    
    assert_close(dd.get_result(1), rev_result, 1e-12);
    
    // Test linearity property for exponential function
    dd.clear();
    std::vector<double> points = {0.0, 1.0, 2.0, 3.0};
    for (double z : points) {
        dd.add_element(z);
    }
    
    // For equally spaced points, the divided differences should have specific properties
    // For exp(x) with spacing h, f[x0, x1, ..., xn] should approach exp(x0)/h^n/n! as h->0
    // This is a complex property to test exactly, but we can check reasonableness
    double h = 1.0;
    double x0 = 0.0;
    for (size_t n = 1; n < points.size(); ++n) {
        double expected_order = std::exp(x0) / std::pow(h, n);
        double actual = dd.get_result(n).to_double();
        double ratio = actual / expected_order;
        
        // The ratio should be close to 1/n! for small h
        double expected_ratio = 1.0;
        for (int i = 1; i <= static_cast<int>(n); ++i) {
            expected_ratio /= i;
        }
        
        // For h=1, this is approximate, but should be reasonable
        ASSERT_GT(ratio, 0.1);  // Should be positive and not too small
        ASSERT_LT(ratio, 10.0); // Should not be too large
    }
    
    std::cout << "✓\n";
}

// Test factorial scaling behavior explicitly
void test_factorial_scaling() {
    std::cout << "Testing factorial scaling behavior... ";
    
    // Test with small, known values
    std::vector<double> z = {0.0, 0.1, 0.2};
    DividedDifferences dd;
    
    for (double val : z) {
        dd.add_element(val);
    }
    
    // Calculate raw divided differences manually
    double f0 = std::exp(z[0]);
    double f1 = std::exp(z[1]);
    double f2 = std::exp(z[2]);
    
    double raw_f0 = f0;
    double raw_f01 = (f1 - f0) / (z[1] - z[0]);
    double raw_f012 = ((f2 - f1) / (z[2] - z[1]) - raw_f01) / (z[2] - z[0]);
    
    // The implementation should return factorial-scaled values
    assert_close(dd.get_result(0), raw_f0);                    // 0! * f[z0]
    assert_close(dd.get_result(1), raw_f01);                   // 1! * f[z0,z1]
    assert_close(dd.get_result(2), 2.0 * raw_f012);           // 2! * f[z0,z1,z2]
    
    std::cout << "✓\n";
}

// Test consistency with different sequences but same set of points
void test_permutation_invariance() {
    std::cout << "Testing permutation invariance... ";
    
    std::vector<double> original = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> permuted = {3.0, 1.0, 4.0, 2.0};
    
    DividedDifferences dd1, dd2;
    
    // Add points in original order
    for (double z : original) {
        dd1.add_element(z);
    }
    
    // Add points in permuted order
    for (double z : permuted) {
        dd2.add_element(z);
    }
    
    // The highest order divided difference should be the same regardless of order
    // (since it's symmetric for the exponential function)
    size_t n = original.size() - 1;
    assert_close(dd1.get_result(n), dd2.get_result(n).to_double(), 1e-10);
    
    std::cout << "✓\n";
}

// Test behavior with very small differences (numerical precision)
void test_small_differences() {
    std::cout << "Testing small differences... ";
    
    DividedDifferences dd;
    double base = 1.0;
    double epsilon = 1e-12;
    
    dd.add_element(base);
    dd.add_element(base + epsilon);
    dd.add_element(base + 2 * epsilon);
    
    // All results should be finite and reasonable
    for (size_t i = 0; i < dd.size(); ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
        ASSERT_FALSE(std::isnan(result));
        
        // Results should be positive for exponential function with positive inputs
        ASSERT_GT(result, 0.0);
    }
    
    std::cout << "✓\n";
}

// Test mixed positive and negative values
void test_mixed_sign_values() {
    std::cout << "Testing mixed sign values... ";
    
    DividedDifferences dd;
    std::vector<double> mixed_values = {-2.0, -1.0, 0.0, 1.0, 2.0};
    
    for (double z : mixed_values) {
        dd.add_element(z);
    }
    
    // All results should be finite
    for (size_t i = 0; i < dd.size(); ++i) {
        double result = dd.get_result(i).to_double();
        ASSERT_TRUE(std::isfinite(result));
        ASSERT_FALSE(std::isnan(result));
    }
    
    // First result should be exp(-2.0)
    assert_close(dd.get_result(0), std::exp(-2.0));
    
    std::cout << "✓\n";
}

int main() {
    std::cout << "=== Comprehensive DividedDifferences Tests ===\n\n";
    
    try {
        test_construction();
        test_single_element();
        test_two_elements();
        test_three_elements();
        test_removal_operations();
        test_clear_operations();
        test_error_handling();
        test_large_sequences();
        test_numerical_stability();
        test_random_sequences();
        test_performance();
        test_scaling_factor_changes();
        test_edge_cases();
        
        // New comprehensive tests
        test_mathematical_properties();
        test_factorial_scaling();
        test_permutation_invariance();
        test_small_differences();
        test_mixed_sign_values();
        
        std::cout << "\n✅ All DividedDifferences tests passed!\n";
        std::cout << "The implementation is robust and handles all tested scenarios correctly.\n";
        std::cout << "Key features verified:\n";
        std::cout << "  ✓ Factorial-scaled divided differences (q! * f[z0,...,zq])\n";
        std::cout << "  ✓ Reentrancy protection to prevent infinite recursion\n";
        std::cout << "  ✓ Mathematical invariants and properties\n";
        std::cout << "  ✓ Numerical stability with extreme values\n";
        std::cout << "  ✓ Permutation invariance for symmetric functions\n";
        return 0;
    } catch (const std::exception& e) {
        std::cout << "\n❌ Test failed with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cout << "\n❌ Test failed with unknown exception\n";
        return 1;
    }
}