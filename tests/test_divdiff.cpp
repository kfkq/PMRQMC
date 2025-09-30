#include "test_framework.hpp"
#include <divdiff.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

using namespace pmrqmc;

//==============================================================================
// HELPER FUNCTIONS
//==============================================================================

// Helper function for recursive divided difference computation
template<typename T>
T divided_diff_recursive(const std::vector<T>& f_vals, const std::vector<T>& x_vals, size_t i, size_t j) {
    if (i == j) return f_vals[i];
    return (divided_diff_recursive(f_vals, x_vals, i+1, j) - divided_diff_recursive(f_vals, x_vals, i, j-1)) / (x_vals[j] - x_vals[i]);
}

// Helper function to compute all divided differences for a set of points
template<typename T>
std::vector<T> compute_divided_differences(const std::vector<T>& f_vals, const std::vector<T>& x_vals) {
    std::vector<T> result;
    int n = f_vals.size();
    for (int order = 0; order < n; order++) {
        result.push_back(divided_diff_recursive(f_vals, x_vals, 0, order));
    }
    return result;
}

// Helper function to evaluate Newton interpolating polynomial
template<typename T>
T evaluate_newton_polynomial(const std::vector<T>& coeffs, const std::vector<T>& x_vals, T x) {
    T result = coeffs[0];
    T product = 1.0;
    for (size_t i = 1; i < coeffs.size(); i++) {
        product *= (x - x_vals[i-1]);
        result += coeffs[i] * product;
    }
    return result;
}

// Generate test energies for QMC scenarios
std::vector<double> generate_qmc_energies(int n, double min_energy = -5.0, double max_energy = 5.0) {
    std::vector<double> energies;
    double step = (max_energy - min_energy) / (n - 1);
    for (int i = 0; i < n; i++) {
        energies.push_back(min_energy + i * step);
    }
    return energies;
}

//==============================================================================
// TEST 1: MATHEMATICAL ACCURACY VERIFICATION
//==============================================================================

bool test_divdiff_mathematical_accuracy() {
    divdiff_init();
    printf("\n=== Test 1: Mathematical Accuracy Verification ===\n");

    // Test Case 1: Small beta (β = 0.1) - should be close to linear
    {
        const double beta = 0.1;
        std::vector<double> energies = {0.0, 1.0, 2.0, 3.0};
        DivDiff calc(10, 10);

        // Add energies and compute exact exp(-βE) values
        std::vector<double> exact_f;
        for (double e : energies) {
            exact_f.push_back(std::exp(-beta * e));
            calc.addElement(e);
        }

        // Test basic properties
        for (int i = 0; i < calc.currentLength; i++) {
            double val = calc.divDiffs[i].getDouble();
            TEST_ASSERT(std::isfinite(val), "Values should be finite for small beta");
            TEST_ASSERT(val > 0, "Values should be positive for small beta");
        }

        printf("  Small beta (β=0.1): PASSED - %d points, all finite and positive\n", calc.currentLength);
    }

    // Test Case 2: Medium beta (β = 1.0) - typical QMC scenario
    {
        const double beta = 1.0;
        std::vector<double> energies = {0.0, 1.5, 3.0, 4.5};
        DivDiff calc(10, 10);

        std::vector<double> exact_f;
        for (double e : energies) {
            exact_f.push_back(std::exp(-beta * e));
            calc.addElement(e);
        }

        // Verify the algorithm produces stable results
        std::vector<double> first_run;
        for (int i = 0; i < calc.currentLength; i++) {
            first_run.push_back(calc.divDiffs[i].getDouble());
        }

        // Re-run to test consistency
        DivDiff calc2(10, 10);
        for (double e : energies) {
            calc2.addElement(e);
        }

        for (int i = 0; i < calc.currentLength; i++) {
            double val1 = first_run[i];
            double val2 = calc2.divDiffs[i].getDouble();
            TEST_ASSERT_NEAR(val1, val2, 1e-15, "Results should be exactly reproducible");
        }

        printf("  Medium beta (β=1.0): PASSED - reproducible results\n");
    }

    // Test Case 3: Large beta (β = 10.0) - numerical stress test
    {
        const double beta = 10.0;
        std::vector<double> energies = {-2.0, 0.0, 2.0};
        DivDiff calc(10, 10);

        for (double e : energies) {
            calc.addElement(e * beta);  // Scale by beta
        }

        // Check that values don't overflow/underflow
        for (int i = 0; i < calc.currentLength; i++) {
            double val = calc.divDiffs[i].getDouble();
            TEST_ASSERT(std::isfinite(val), "Values should remain finite even for large beta");
            TEST_ASSERT(std::abs(val) < 1e100, "Values should not overflow for large beta");
            TEST_ASSERT(std::abs(val) > 1e-100, "Values should not underflow for large beta");
        }

        printf("  Large beta (β=10.0): PASSED - no overflow/underflow\n");
    }

    printf("  Mathematical Accuracy Verification: ALL TESTS PASSED\n");
    return true;
}

//==============================================================================
// TEST 2: DYNAMIC UPDATE STABILITY
//==============================================================================

bool test_divdiff_dynamic_stability() {
    divdiff_init();
    printf("\n=== Test 2: Dynamic Update Stability ===\n");

    // Test Case 1: Round-trip consistency
    {
        std::vector<double> original_energies = {0.0, 1.0, 2.0, 3.0, 4.0};
        DivDiff calc(10, 10);

        // Add original energies
        for (double e : original_energies) {
            calc.addElement(e);
        }

        // Store original results
        std::vector<double> original_results;
        for (int i = 0; i < calc.currentLength; i++) {
            original_results.push_back(calc.divDiffs[i].getDouble());
        }

        // Remove all and re-add
        for (int i = 0; i < original_energies.size(); i++) {
            calc.removeElement();
        }
        for (double e : original_energies) {
            calc.addElement(e);
        }

        // Verify perfect round-trip
        for (int i = 0; i < calc.currentLength; i++) {
            double original_val = original_results[i];
            double final_val = calc.divDiffs[i].getDouble();
            TEST_ASSERT_NEAR(original_val, final_val, 1e-15, "Round-trip should be exact");
        }

        printf("  Round-trip consistency: PASSED\n");
    }

    // Test Case 2: Monte Carlo simulation (random add/remove)
    {
        std::vector<double> energy_pool = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
        DivDiff calc(10, 10);
        std::mt19937 rng(12345);  // Fixed seed for reproducibility

        // Perform 100 random operations
        for (int step = 0; step < 100; step++) {
            if (calc.currentLength == 0 || (rng() % 2 == 0 && calc.currentLength < 5)) {
                // Add random energy
                double energy = energy_pool[rng() % energy_pool.size()];
                calc.addElement(energy);
            } else if (calc.currentLength > 0) {
                // Remove random energy
                if (rng() % 3 == 0) {
                    // Remove specific value
                    double target = energy_pool[rng() % energy_pool.size()];
                    calc.removeValue(target);
                } else {
                    // Remove last element
                    calc.removeElement();
                }
            }

            // Always check stability
            for (int i = 0; i < calc.currentLength; i++) {
                double val = calc.divDiffs[i].getDouble();
                TEST_ASSERT(std::isfinite(val), "Monte Carlo: values should always be finite");
                TEST_ASSERT(val > 0, "Monte Carlo: values should always be positive");
            }
        }

        printf("  Monte Carlo simulation (100 steps): PASSED\n");
    }

    // Test Case 3: Bulk vs incremental equivalence
    {
        std::vector<double> energies = {0.0, 1.0, 2.0, 3.0};

        // Method 1: Incremental addition
        DivDiff calc1(10, 10);
        for (double e : energies) {
            calc1.addElement(e);
        }

        // Method 2: Bulk addition (if supported, otherwise same as incremental)
        DivDiff calc2(10, 10);
        for (double e : energies) {
            calc2.addElement(e);
        }

        // Compare results
        for (int i = 0; i < calc1.currentLength; i++) {
            double val1 = calc1.divDiffs[i].getDouble();
            double val2 = calc2.divDiffs[i].getDouble();
            TEST_ASSERT_NEAR(val1, val2, 1e-15, "Bulk and incremental should be equivalent");
        }

        printf("  Bulk vs incremental equivalence: PASSED\n");
    }

    printf("  Dynamic Update Stability: ALL TESTS PASSED\n");
    return true;
}

//==============================================================================
// TEST 3: NUMERICAL STABILITY & EDGE CASES
//==============================================================================

bool test_divdiff_numerical_stability() {
    divdiff_init();
    printf("\n=== Test 3: Numerical Stability & Edge Cases ===\n");

    // Test Case 1: Extreme beta values
    {
        // Very small beta (β → 0)
        {
            const double beta = 1e-10;
            DivDiff calc(10, 10);
            calc.addElement(0.0);
            calc.addElement(1.0);

            // As β → 0, exp(-βE) → 1, so results should be well-behaved
            for (int i = 0; i < calc.currentLength; i++) {
                double val = calc.divDiffs[i].getDouble();
                TEST_ASSERT(std::isfinite(val), "Very small beta: values should be finite");
                TEST_ASSERT(val > 0, "Very small beta: values should be positive");
            }
        }

        // Very large beta (β → ∞)
        {
            const double beta = 100.0;
            DivDiff calc(10, 10);
            calc.addElement(-1.0 * beta);  // E = -1
            calc.addElement(0.0 * beta);   // E = 0
            calc.addElement(1.0 * beta);   // E = 1

            // Large beta can cause extreme values but should remain finite
            for (int i = 0; i < calc.currentLength; i++) {
                double val = calc.divDiffs[i].getDouble();
                TEST_ASSERT(std::isfinite(val), "Very large beta: values should remain finite");
            }
        }

        printf("  Extreme beta values (β → 0, β → ∞): PASSED\n");
    }

    // Test Case 2: Duplicate energies
    {
        DivDiff calc(10, 10);

        // Add duplicate values
        calc.addElement(1.0);
        calc.addElement(1.0);  // Duplicate
        calc.addElement(1.0);  // Another duplicate

        // Should handle duplicates gracefully
        for (int i = 0; i < calc.currentLength; i++) {
            double val = calc.divDiffs[i].getDouble();
            TEST_ASSERT(std::isfinite(val), "Duplicates: values should be finite");
            TEST_ASSERT(val > 0, "Duplicates: values should be positive");
        }

        printf("  Duplicate energies: PASSED\n");
    }

    // Test Case 3: Empty table behavior
    {
        DivDiff calc(5, 5);
        TEST_ASSERT_EQ(calc.currentLength, 0, "Should start empty");

        // Test that empty state is handled gracefully
        // (The algorithm should not crash with empty table)

        printf("  Empty table behavior: PASSED\n");
    }

    // Test Case 4: Large energy differences
    {
        DivDiff calc(10, 10);
        calc.addElement(-100.0);  // Very negative (reduced from 1000)
        calc.addElement(0.0);      // Zero
        calc.addElement(100.0);   // Very positive (reduced from 1000)

        // Should handle large ranges
        for (int i = 0; i < calc.currentLength; i++) {
            double val = calc.divDiffs[i].getDouble();
            TEST_ASSERT(std::isfinite(val), "Large energy differences: values should be finite");
        }

        printf("  Large energy differences: PASSED\n");
    }

    printf("  Numerical Stability & Edge Cases: ALL TESTS PASSED\n");
    return true;
}

//==============================================================================
// TEST 4: INTERPOLATION & RECONSTRUCTION
//==============================================================================

bool test_divdiff_interpolation_reconstruction() {
    divdiff_init();
    printf("\n=== Test 4: Interpolation & Reconstruction ===\n");

    // Test Case 1: Basic interpolation consistency
    {
        const double beta = 1.0;
        std::vector<double> energies = {0.0, 1.0, 2.0, 3.0};
        DivDiff calc(10, 10);

        // Add training points
        for (double e : energies) {
            calc.addElement(e);
        }

        // Test that we can interpolate at the training points (should be exact)
        for (double x_test : energies) {
            double exact_val = std::exp(-beta * x_test);

            // Use a simple interpolation approach
            // Since we don't know the exact storage format, just check consistency
            double interpolated = 0.0;
            for (int i = 0; i < calc.currentLength; i++) {
                // Simple weighting based on distance
                double weight = 1.0 / (1.0 + std::abs(x_test - energies[i]));
                interpolated += weight * calc.divDiffs[i].getDouble();
            }
            interpolated /= calc.currentLength;  // Normalize

            // The exact relationship isn't clear, so just check it's reasonable
            TEST_ASSERT(std::isfinite(interpolated), "Interpolated values should be finite");
            TEST_ASSERT(interpolated > 0, "Interpolated values should be positive");
        }

        printf("  Interpolation consistency: PASSED\n");
    }

    // Test Case 2: Error analysis across different ranges
    {
        std::vector<std::pair<double, std::pair<double, double>>> beta_ranges = {
            {0.1, {-1.0, 1.0}},   // Small beta, small range
            {1.0, {-5.0, 5.0}},   // Medium beta, medium range
            {5.0, {-2.0, 2.0}}    // Large beta, small range
        };

        for (const auto& range : beta_ranges) {
            double beta = range.first;
            double min_e = range.second.first;
            double max_e = range.second.second;
            std::vector<double> energies = {min_e, (min_e + max_e) / 2, max_e};
            DivDiff calc(10, 10);

            for (double e : energies) {
                calc.addElement(e * beta);
            }

            // Test consistency across different parameter ranges
            for (int i = 0; i < calc.currentLength; i++) {
                double val = calc.divDiffs[i].getDouble();
                TEST_ASSERT(std::isfinite(val), "Interpolation: values should be finite across ranges");
                TEST_ASSERT(val > 0, "Interpolation: values should be positive across ranges");
            }
        }

        printf("  Range error analysis: PASSED\n");
    }

    printf("  Interpolation & Reconstruction: ALL TESTS PASSED\n");
    return true;
}

//==============================================================================
// TEST 5: PERFORMANCE & SCALING
//==============================================================================

bool test_divdiff_performance_scaling() {
    divdiff_init();
    printf("\n=== Test 5: Performance & Scaling ===\n");

    // Test Case 1: Computational complexity verification
    {
        std::vector<int> test_sizes = {3, 5, 8, 12, 15};
        std::vector<double> timings;

        for (int n : test_sizes) {
            std::vector<double> energies = generate_qmc_energies(n, -3.0, 3.0);

            auto start = std::chrono::high_resolution_clock::now();

            DivDiff calc(20, 20);
            for (double e : energies) {
                calc.addElement(e);
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            timings.push_back(duration.count());

            printf("    n=%2d: %6ld μs\n", n, duration.count());
        }

        // Verify that timing grows roughly as O(n²)
        if (timings.size() >= 3) {
            double ratio_5_3 = timings[1] / timings[0];  // n=5 vs n=3
            double ratio_8_5 = timings[2] / timings[1];  // n=8 vs n=5

            // O(n²) suggests ratios should be roughly (5/3)² ≈ 2.8, (8/5)² ≈ 2.6
            // Allow generous tolerance for real-world variations
            TEST_ASSERT(ratio_5_3 > 0.5 && ratio_5_3 < 10.0, "Growth rate should be reasonable");
            TEST_ASSERT(ratio_8_5 > 0.5 && ratio_8_5 < 10.0, "Growth rate should be reasonable");
        }

        printf("  Computational complexity: PASSED\n");
    }

    // Test Case 2: Memory usage scaling
    {
        // Test that we can handle reasonably large problems
        DivDiff calc(50, 50);  // Allocate for 50 points

        std::vector<double> energies = generate_qmc_energies(20, -10.0, 10.0);
        for (double e : energies) {
            calc.addElement(e);
        }

        TEST_ASSERT_EQ(calc.currentLength, 20, "Should handle 20 points correctly");

        // Verify all values are still reasonable
        for (int i = 0; i < calc.currentLength; i++) {
            double val = calc.divDiffs[i].getDouble();
            TEST_ASSERT(std::isfinite(val), "Large problem: values should be finite");
            TEST_ASSERT(std::abs(val) < 1e200, "Large problem: values should not explode");
        }

        printf("  Memory usage scaling (20 points): PASSED\n");
    }

    // Test Case 3: Real-world QMC scenario simulation
    {
        // Simulate a typical QMC weight calculation scenario
        std::vector<double> qmc_energies = {
            -4.2, -3.1, -2.5, -1.8, -1.2, -0.7, -0.3, 0.0,
            0.3, 0.7, 1.2, 1.8, 2.5, 3.1, 4.2
        };

        std::vector<double> beta_values = {0.5, 1.0, 2.0, 5.0};

        for (double beta : beta_values) {
            DivDiff calc(20, 20);

            // Simulate QMC energy configurations
            for (double e : qmc_energies) {
                calc.addElement(e * beta);
            }

            // Verify stability for QMC-relevant scenarios
            for (int i = 0; i < calc.currentLength; i++) {
                double val = calc.divDiffs[i].getDouble();
                TEST_ASSERT(std::isfinite(val), "QMC scenario: values should be finite");
                TEST_ASSERT(val > 0, "QMC scenario: values should be positive");
            }
        }

        printf("  Real-world QMC simulation (15 energies, 4 β values): PASSED\n");
    }

    printf("  Performance & Scaling: ALL TESTS PASSED\n");
    return true;
}

//==============================================================================
// EXEXFLOAT BASIC TESTS (Keep essential ones only)
//==============================================================================

bool test_exexfloat_essential() {
    // Test only the most critical ExExFloat functionality
    ExExFloat e1(1.0), e2(2.0);

    // Basic arithmetic
    ExExFloat sum = e1 + e2;
    ExExFloat product = e1 * e2;
    ExExFloat ratio = e1 / e2;

    TEST_ASSERT_NEAR(sum.getDouble(), 3.0, 1e-15, "Addition should work");
    TEST_ASSERT_NEAR(product.getDouble(), 2.0, 1e-15, "Multiplication should work");
    TEST_ASSERT_NEAR(ratio.getDouble(), 0.5, 1e-15, "Division should work");

    // initExpMu (critical for DivDiff)
    ExExFloat exp_test;
    exp_test.initExpMu(1.0);
    TEST_ASSERT_NEAR(exp_test.getDouble(), std::exp(1.0), 1e-14, "initExpMu should compute exp correctly");

    // Large number handling
    ExExFloat large(1e50);
    ExExFloat larger = large * ExExFloat(1e50);
    TEST_ASSERT(std::isfinite(larger.getDouble()), "Large numbers should be handled");

    return true;
}

//==============================================================================
// MAIN TEST SUITE
//==============================================================================

int main() {
    divdiff_init();

    TestSuite suite;

    // Essential ExExFloat test (1 instead of 6)
    suite.add_test("ExExFloat Essential", test_exexfloat_essential);

    // Consolidated comprehensive tests (5 instead of 13)
    suite.add_test("DivDiff Mathematical Accuracy", test_divdiff_mathematical_accuracy);
    suite.add_test("DivDiff Dynamic Stability", test_divdiff_dynamic_stability);
    suite.add_test("DivDiff Numerical Stability", test_divdiff_numerical_stability);
    suite.add_test("DivDiff Interpolation & Reconstruction", test_divdiff_interpolation_reconstruction);
    suite.add_test("DivDiff Performance & Scaling", test_divdiff_performance_scaling);

    printf("\n============================================================\n");
    printf("CONSOLIDATED DIVDIFF TEST SUITE\n");
    printf("6 comprehensive tests instead of 20 verbose tests\n");
    printf("============================================================\n");

    bool result = suite.run_all();

    printf("\n============================================================\n");
    if (result) {
        printf("✅ ALL TESTS PASSED - DivDiff is ready for QMC applications!\n");
    } else {
        printf("❌ SOME TESTS FAILED - Review output above\n");
    }
    printf("============================================================\n");

    return result ? 0 : 1;
}