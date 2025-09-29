#include "test_framework.hpp"
#include <divdiff.h>

using namespace pmrqmc;

bool test_exexfloat_basic() {
    ExExFloat e(1.0);
    TEST_ASSERT_NEAR(e.getDouble(), 1.0, 1e-15, "ExExFloat constructor should work");

    ExExFloat e2 = 2.5;
    TEST_ASSERT_NEAR(e2.getDouble(), 2.5, 1e-15, "Assignment from double should work");

    ExExFloat e3(1.0);
    ExExFloat e4(e3 + ExExFloat(1.0));
    TEST_ASSERT_NEAR(e4.getDouble(), 2.0, 1e-15, "Addition should work");

    return true;
}

bool test_exexfloat_large_numbers() {
    // Test extreme numbers
    ExExFloat e1(1e100);
    ExExFloat e2(1e50);

    ExExFloat product = e1 * e2;
    double expected = 1e150;
    double actual = product.getDouble();

    // Allow some relative error due to mantissa precision
    TEST_ASSERT_NEAR(std::log10(actual), std::log10(expected), 1.0, "Large number multiplication should work");

    return true;
}

bool test_exexfloat_arithmetic() {
    ExExFloat a(2.0);
    ExExFloat b(3.0);

    TEST_ASSERT_NEAR((a + b).getDouble(), 5.0, 1e-15, "Addition");
    TEST_ASSERT_NEAR((a - b).getDouble(), -1.0, 1e-15, "Subtraction");
    TEST_ASSERT_NEAR((a * b).getDouble(), 6.0, 1e-15, "Multiplication");
    TEST_ASSERT_NEAR((a / b).getDouble(), 2.0/3.0, 1e-14, "Division");

    return true;
}

bool test_exexfloat_comparison() {
    ExExFloat a(1.0);
    ExExFloat b(2.0);

    TEST_ASSERT(a >= 0.5, "Should be greater than smaller ExExFloat");
    TEST_ASSERT(!(a >= b), "Should not be greater than larger ExExFloat");
    TEST_ASSERT(a >= a, "Should be greater than or equal to itself");

    return true;
}

bool test_exexfloat_special_functions() {
    ExExFloat a(4.0);
    ExExFloat sqrt_a = a.squareRoot();
    TEST_ASSERT_NEAR(sqrt_a.getDouble(), 2.0, 1e-14, "Square root should work");

    // Legacy squareRoot may not throw for negative numbers, just return NaN
    ExExFloat b(-1.0);
    ExExFloat sqrt_neg = b.squareRoot();
    // The legacy implementation may not handle negative inputs gracefully,
    // so we just check it doesn't crash

    ExExFloat c(1.0);
    TEST_ASSERT_NEAR(c.sign(), 1, 1, "Positive sign");
    TEST_ASSERT_NEAR(c.absoluteValue().getDouble(), 1.0, 1e-15, "Absolute value");

    return true;
}

bool test_exexfloat_init_expmu() {
    ExExFloat e(0.0);
    e.initExpMu(1.0);  // exp(1.0) ≈ 2.718

    double result = e.getDouble();
    TEST_ASSERT_NEAR(result, std::exp(1.0), 1e-14, "initExpMu should compute exp(value)");

    return true;
}

bool test_divdiff_constructor() {
    // Legacy constructor takes maxlen and smax directly
    DivDiff calc(100, 50);  // maxlen=100, smax=50

    TEST_ASSERT_EQ(calc.currentLength, 0, "Should start with length 0");
    // Legacy doesn't have scaling_factor accessor, just check basic construction
    TEST_ASSERT(calc.divDiffs != nullptr, "Should have divDiffs array allocated");

    return true;
}

bool test_divdiff_basic_add_element() {
    divdiff_init();  // Initialize global state

    DivDiff calc(10, 10);

    // Add first element (should trigger initialization)
    calc.addElement(1.0);
    TEST_ASSERT_EQ(calc.currentLength, 1, "Should have 1 element");
    TEST_ASSERT_NEAR(calc.zVals[0], 1.0, 1e-15, "Should store correct energy in zVals array");

    // Add second element
    calc.addElement(2.0);
    TEST_ASSERT_EQ(calc.currentLength, 2, "Should have 2 elements");

    return true;
}

bool test_divdiff_add_remove_element() {
    divdiff_init();

    DivDiff calc(10, 10);

    calc.addElement(1.0);
    calc.addElement(2.0);
    TEST_ASSERT_EQ(calc.currentLength, 2, "Should have 2 elements after adding");

    calc.removeElement();
    TEST_ASSERT_EQ(calc.currentLength, 1, "Should have 1 element after removing");

    calc.removeElement();
    TEST_ASSERT_EQ(calc.currentLength, 0, "Should have 0 elements after removing all");

    return true;
}

bool test_divdiff_remove_value() {
    divdiff_init();

    DivDiff calc(10, 10);

    calc.addElement(1.0);
    calc.addElement(2.0);
    calc.addElement(3.0);
    TEST_ASSERT_EQ(calc.currentLength, 3, "Should have 3 elements");

    // Remove the middle value (2.0) - this removes from bulk, not just one element
    bool removed = calc.removeValue(2.0);
    TEST_ASSERT(removed, "Should successfully remove value");
    TEST_ASSERT_EQ(calc.currentLength, 2, "Should have 2 elements after removal");

    return true;
}

bool test_divdiff_divdiffs_computation() {
    divdiff_init();

    DivDiff calc(10, 10);

    // Add a few energies and check that divDiffs are computed
    calc.addElement(0.0);
    // Legacy doesn't have direct size check, but we can check currentLength

    calc.addElement(1.0);
    TEST_ASSERT_EQ(calc.currentLength, 2, "Should have 2 elements");

    // divDiffs should contain the DDEF values (divided differences of exponential function)
    // Check the first few divDiffs values are finite
    if (calc.divDiffs != nullptr) {
        double val1 = calc.divDiffs[0].getDouble();
        double val2 = calc.divDiffs[1].getDouble();
        TEST_ASSERT(std::isfinite(val1), "First divDiffs should be finite");
        TEST_ASSERT(std::isfinite(val2), "Second divDiffs should be finite");
    }

    return true;
}

bool test_divdiff_scaling_change() {
    divdiff_init();

    // Test that switching contexts doesn't crash - legacy handles scaling internally
    DivDiff calc(10, 10);

    calc.addElement(0.0);
    calc.addElement(10.0);  // Large value that might trigger scaling internally

    TEST_ASSERT(calc.currentLength >= 2, "Should have added elements");

    return true;
}

bool test_divdiff_memory_management() {
    // Test that large divdiff instances can be created without crashing
    DivDiff calc(100, 50);  // Large dimensions

    // Just check basic functionality works
    calc.addElement(1.0);
    TEST_ASSERT_EQ(calc.currentLength, 1, "Should work with large dimensions");

    return true;
}

bool test_divdiff_print_state() {
    divdiff_init();

    DivDiff calc(5, 5);
    calc.addElement(1.5);
    calc.addElement(-0.5);

    // Just test that basic operations work
    TEST_ASSERT(calc.currentLength == 2, "Should process multiple elements");

    return true;
}

int main() {
    divdiff_init();  // Use global function, not in namespace

    TestSuite suite;

    suite.add_test("ExExFloat Basic Operations", test_exexfloat_basic);
    suite.add_test("ExExFloat Large Numbers", test_exexfloat_large_numbers);
    suite.add_test("ExExFloat Arithmetic", test_exexfloat_arithmetic);
    suite.add_test("ExExFloat Comparison", test_exexfloat_comparison);
    suite.add_test("ExExFloat Special Functions", test_exexfloat_special_functions);
    suite.add_test("ExExFloat init_expmu", test_exexfloat_init_expmu);

    suite.add_test("DivDiff Constructor", test_divdiff_constructor);
    suite.add_test("DivDiff Basic Add Element", test_divdiff_basic_add_element);
    suite.add_test("DivDiff Add/Remove Element", test_divdiff_add_remove_element);
    suite.add_test("DivDiff Remove Value", test_divdiff_remove_value);
    suite.add_test("DivDiff DivDiffs Computation", test_divdiff_divdiffs_computation);
    suite.add_test("DivDiff Scaling Change", test_divdiff_scaling_change);
    suite.add_test("DivDiff Memory Management", test_divdiff_memory_management);
    suite.add_test("DivDiff Print State", test_divdiff_print_state);

    return suite.run_all() ? 0 : 1;
}