#include <iostream>
#include <cmath>
#include <cassert>
#include "core/ExExFloat.hpp"

using namespace pmrqmc::core;

void test_basic_construction() {
    std::cout << "Testing basic construction... ";
    
    ExExFloat a;
    assert(a.to_double() == 0.0);
    assert(a.sgn() == 0);
    
    ExExFloat b(1.5);
    assert(std::abs(b.to_double() - 1.5) < 1e-15);
    assert(b.sgn() == 1);
    
    ExExFloat c(-2.0);
    assert(std::abs(c.to_double() - (-2.0)) < 1e-15);
    assert(c.sgn() == -1);
    
    std::cout << "PASS\n";
}

void test_arithmetic() {
    std::cout << "Testing arithmetic operations... ";
    
    ExExFloat a(1.0);
    ExExFloat b(2.0);
    ExExFloat c = a + b;
    assert(std::abs(c.to_double() - 3.0) < 1e-15);
    
    ExExFloat d(5.0);
    ExExFloat e(3.0);
    ExExFloat f = d - e;
    assert(std::abs(f.to_double() - 2.0) < 1e-15);
    
    ExExFloat g(3.0);
    ExExFloat h(4.0);
    ExExFloat i = g * h;
    assert(std::abs(i.to_double() - 12.0) < 1e-15);
    
    ExExFloat j(10.0);
    ExExFloat k(2.0);
    ExExFloat l = j / k;
    assert(std::abs(l.to_double() - 5.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_mixed_arithmetic() {
    std::cout << "Testing mixed arithmetic with doubles... ";
    
    ExExFloat a(2.0);
    
    assert(std::abs((a + 3.0).to_double() - 5.0) < 1e-15);
    assert(std::abs((5.0 + a).to_double() - 7.0) < 1e-15);
    assert(std::abs((a - 1.0).to_double() - 1.0) < 1e-15);
    assert(std::abs((8.0 - a).to_double() - 6.0) < 1e-15);
    assert(std::abs((a * 3.0).to_double() - 6.0) < 1e-15);
    assert(std::abs((4.0 * a).to_double() - 8.0) < 1e-15);
    assert(std::abs((a / 2.0).to_double() - 1.0) < 1e-15);
    assert(std::abs((12.0 / a).to_double() - 6.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_compound_assignment() {
    std::cout << "Testing compound assignment... ";
    
    ExExFloat a(1.0);
    ExExFloat b(2.0);
    a += b;
    assert(std::abs(a.to_double() - 3.0) < 1e-15);
    
    ExExFloat c(5.0);
    ExExFloat d(3.0);
    c -= d;
    assert(std::abs(c.to_double() - 2.0) < 1e-15);
    
    ExExFloat e(3.0);
    ExExFloat f(4.0);
    e *= f;
    assert(std::abs(e.to_double() - 12.0) < 1e-15);
    
    ExExFloat g(10.0);
    ExExFloat h(2.0);
    g /= h;
    assert(std::abs(g.to_double() - 5.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_mathematical_functions() {
    std::cout << "Testing mathematical functions... ";
    
    ExExFloat a(3.0);
    ExExFloat b(-3.0);
    assert(std::abs(a.abs().to_double() - 3.0) < 1e-15);
    assert(std::abs(b.abs().to_double() - 3.0) < 1e-15);
    
    ExExFloat c(4.0);
    ExExFloat d = c.sqrt();
    assert(std::abs(d.to_double() - 2.0) < 1e-15);
    
    ExExFloat e(0.0);
    ExExFloat f = e.sqrt();
    assert(f.to_double() == 0.0);
    
    ExExFloat g(3.0);
    ExExFloat h = g.square();
    assert(std::abs(h.to_double() - 9.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_from_exp() {
    std::cout << "Testing from_exp factory function... ";
    
    ExExFloat a = ExExFloat::from_exp(0.0);  // exp(0) = 1
    assert(std::abs(a.to_double() - 1.0) < 1e-15);
    
    ExExFloat b = ExExFloat::from_exp(std::log(2.0));  // exp(ln(2)) = 2
    assert(std::abs(b.to_double() - 2.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_comparisons() {
    std::cout << "Testing comparison operators... ";
    
    ExExFloat a(3.0);
    ExExFloat b(3.0);
    ExExFloat c(5.0);
    ExExFloat d(-1.0);
    
    assert(a == b);
    assert(a != c);
    assert(a < c);
    assert(c > a);
    assert(a > d);
    
    std::cout << "PASS\n";
}

void test_utility_functions() {
    std::cout << "Testing utility functions... ";
    
    ExExFloat a(0.0);
    assert(a.is_zero());
    assert(a.is_finite());
    assert(!a.is_nan());
    assert(!a.is_inf());
    
    ExExFloat b(3.0);
    assert(!b.is_zero());
    assert(b.is_finite());
    assert(!b.is_nan());
    assert(!b.is_inf());
    
    ExExFloat c(std::numeric_limits<double>::infinity());
    assert(!c.is_zero());
    assert(!c.is_finite());
    assert(!c.is_nan());
    assert(c.is_inf());
    
    ExExFloat d(std::numeric_limits<double>::quiet_NaN());
    assert(!d.is_zero());
    assert(!d.is_finite());
    assert(d.is_nan());
    assert(!d.is_inf());
    
    std::cout << "PASS\n";
}

void test_cache_performance() {
    std::cout << "Testing cache functionality... ";
    
    ExExFloat a(1.0);
    ExExFloat b(2.0);
    
    // Multiple operations that should benefit from cache
    for (int i = 0; i < 10; ++i) {
        a += b;
    }
    
    // Result should be 1 + 10*2 = 21
    assert(std::abs(a.to_double() - 21.0) < 1e-15);
    
    std::cout << "PASS\n";
}

void test_extreme_values() {
    std::cout << "Testing extreme values... ";
    
    // Test very large numbers
    ExExFloat a(1e300);
    ExExFloat b(1e300);
    ExExFloat c = a + b;
    assert(c.to_double() > 1e300);
    
    // Test very small numbers
    ExExFloat d(1e-300);
    ExExFloat e(1e-300);
    ExExFloat f = d + e;
    assert(f.to_double() > 0.0);
    
    // Test zero handling
    ExExFloat zero(0.0);
    ExExFloat num(5.0);
    assert(std::abs((zero + num).to_double() - 5.0) < 1e-15);
    assert(std::abs((num + zero).to_double() - 5.0) < 1e-15);
    
    std::cout << "PASS\n";
}

int main() {
    std::cout << "Running ExExFloat tests...\n\n";
    
    test_basic_construction();
    test_arithmetic();
    test_mixed_arithmetic();
    test_compound_assignment();
    test_mathematical_functions();
    test_from_exp();
    test_comparisons();
    test_utility_functions();
    test_cache_performance();
    test_extreme_values();
    
    std::cout << "\nAll tests passed! ExExFloat is working correctly.\n";
    return 0;
}