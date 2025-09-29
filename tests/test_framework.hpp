#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <functional>

namespace pmrqmc {

// Simple test framework
class TestSuite {
public:
    struct TestResult {
        std::string name;
        bool passed;
        std::string message;
    };

    void add_test(const std::string& name, std::function<bool()> test_func) {
        tests_.emplace_back(name, test_func);
    }

    bool run_all() {
        int passed = 0;
        int failed = 0;
        std::vector<TestResult> results;

        std::cout << "Running " << tests_.size() << " tests..." << std::endl;
        std::cout << std::string(50, '=') << std::endl;

        for (const auto& [name, test_func] : tests_) {
            TestResult result{name, false, ""};

            try {
                result.passed = test_func();
                if (result.passed) {
                    passed++;
                } else {
                    failed++;
                    result.message = "Test function returned false";
                }
            } catch (const std::exception& e) {
                failed++;
                result.message = std::string("Exception: ") + e.what();
            } catch (...) {
                failed++;
                result.message = "Unknown exception";
            }

            results.push_back(result);

            if (result.passed) {
                std::cout << "[PASS] ";
            } else {
                std::cout << "[FAIL] ";
            }
            std::cout << name << std::endl;

            if (!result.passed && !result.message.empty()) {
                std::cout << "       Reason: " << result.message << std::endl;
            }
        }

        std::cout << std::string(50, '=') << std::endl;
        std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;

        if (failed > 0) {
            std::cout << "\nFailed tests:" << std::endl;
            for (const auto& result : results) {
                if (!result.passed) {
                    std::cout << "- " << result.name;
                    if (!result.message.empty()) {
                        std::cout << ": " << result.message;
                    }
                    std::cout << std::endl;
                }
            }
        }

        return failed == 0;
    }

private:
    std::vector<std::pair<std::string, std::function<bool()>>> tests_;
};

// Helper macros for assertions
#define TEST_ASSERT(expr, message) \
    if (!(expr)) { \
        throw std::runtime_error(message); \
    }

#define TEST_ASSERT_EQ(a, b, message) \
    if ((a) != (b)) { \
        throw std::runtime_error(message); \
    }

#define TEST_ASSERT_NEAR(a, b, tol, message) \
    if (std::abs((a) - (b)) > (tol)) { \
        throw std::runtime_error(message); \
    }

} // namespace pmrqmc