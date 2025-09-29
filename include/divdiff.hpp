#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace pmrqmc {

/**
 * Extended-exponent floating-point number for handling very large ranges
 * in divided differences of exponential function computations.
 *
 * Uses a double-precision mantissa with extended integer exponent range.
 */
class ExExFloat {
public:
    // Maximum exponent value for safety and precomputed powers
    static constexpr int max_exponent = 100'000;

    // Constructors
    ExExFloat() noexcept = default;
    ExExFloat(double val) noexcept;
    ExExFloat(const ExExFloat&) noexcept = default;
    ExExFloat(ExExFloat&&) noexcept = default;

    // Assignment
    ExExFloat& operator=(const ExExFloat&) noexcept = default;
    ExExFloat& operator=(ExExFloat&&) noexcept = default;
    ExExFloat& operator=(double rhs) noexcept;

    // Arithmetic operators
    ExExFloat operator+(const ExExFloat& rhs) const noexcept;
    ExExFloat operator-(const ExExFloat& rhs) const noexcept;
    ExExFloat operator*(const ExExFloat& rhs) const noexcept;
    ExExFloat operator/(const ExExFloat& rhs) const;

    ExExFloat operator*(double rhs) const noexcept;
    ExExFloat operator/(double rhs) const;
    friend ExExFloat operator*(double lhs, const ExExFloat& rhs) noexcept;
    friend ExExFloat operator/(double lhs, const ExExFloat& rhs);

    // Compound assignment
    ExExFloat& operator+=(const ExExFloat& rhs) noexcept;
    ExExFloat& operator-=(const ExExFloat& rhs) noexcept;
    ExExFloat& operator*=(const ExExFloat& rhs) noexcept;
    ExExFloat& operator/=(const ExExFloat& rhs);
    ExExFloat& operator*=(double rhs) noexcept;
    ExExFloat& operator/=(double rhs);

    // Comparison (assumes non-negative mantissas)
    bool operator>=(const ExExFloat& rhs) const noexcept;
    bool operator>=(double rhs) const noexcept;

    // Utility functions
    [[nodiscard]] double double_value() const noexcept;
    [[nodiscard]] int sign() const noexcept;
    [[nodiscard]] ExExFloat abs() const noexcept;
    [[nodiscard]] ExExFloat sqrt() const;

    void init_expmu(double mu);
    void print(std::ostream& os = std::cout) const;

private:
    double mantissa{0.5};
    int exponent{0};

    void normalize() noexcept;

    friend std::ostream& operator<<(std::ostream& os, const ExExFloat& val) {
        val.print(os);
        return os;
    }
};

/**
 * Efficient calculator for divided differences of the exponential function.
 * Used for computing QMC configuration weights.
 */
class DivDiff {
public:
    // Configuration builder
    class Builder {
    public:
        Builder& max_length(int n);
        Builder& max_scaling_factor(int s);

        [[nodiscard]] DivDiff build() const;
        [[nodiscard]] std::unique_ptr<DivDiff> build_unique() const noexcept;

    private:
        int max_len_ = 10'001;
        int max_scaling_ = 500;
    };

    // Constructors
    explicit DivDiff(int max_length, int max_scaling_factor = 500);
    DivDiff(const DivDiff&) = delete;
    DivDiff(DivDiff&&) noexcept = default;
    DivDiff& operator=(const DivDiff&) = delete;
    DivDiff& operator=(DivDiff&&) noexcept = default;
    ~DivDiff() = default;

    // Core interface methods
    void add_element(double energy);
    void remove_element();
    bool remove_value(double value);  // Remove specific value from middle

    // Accessors
    [[nodiscard]] const std::vector<ExExFloat>& divdiffs() const noexcept;
    [[nodiscard]] const std::vector<double>& energies() const noexcept;
    [[nodiscard]] int current_length() const noexcept;
    [[nodiscard]] int scaling_factor() const noexcept;

    // Utility
    [[nodiscard]] size_t memory_usage() const noexcept;
    void print_state(std::ostream& os = std::cout) const;

private:
    // Data members
    int max_length_;
    int max_scaling_;
    int current_length_{0};
    int scaling_factor_{1};
    double mean_energy_{0.0};

    std::vector<double> energies_;
    std::vector<ExExFloat> h_matrix_;
    std::vector<ExExFloat> divdiffs_;
    std::vector<ExExFloat> ddd_matrix_;  // Flattened 2D matrix

    // Helper methods
    void initialize_single_point();
    void add_all_elements(int len);
    [[nodiscard]] bool needs_scaling_change() const noexcept;
    void compute_new_element();
    void update_matrices();

    // Memory management
    void allocate_memory();
    void free_memory() noexcept;
    void resize_if_needed();
};

// Global functions
void divdiff_init();
void divdiff_cleanup() noexcept;

} // namespace pmrqmc