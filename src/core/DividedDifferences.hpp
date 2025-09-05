#pragma once

#include "ExExFloat.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pmrqmc::core {

/**
 * @class DividedDifferences
 * @brief Computes divided differences of the exponential function (DDEF) efficiently.
 *
 * This class implements the algorithm described in:
 * L. Gupta, L. Barash, I. Hen, "Calculating the divided differences of the
 * exponential function by addition and removal of inputs",
 * Computer Physics Communications 254, 107385 (2020).
 *
 * It allows for the sequential addition and removal of input points (z_i)
 * while updating the DDEF vector exp[z_0, ..., z_k] for all k up to the
 * current size. The implementation is memory-safe, using RAII principles
 * and std::vector, and is designed for high performance in QMC simulations.
 */
class DividedDifferences {
public:
    /**
     * @brief Constructs a DividedDifferences calculator.
     * @param max_len Initial capacity for the number of input points.
     * @param s_max Initial capacity for the scaling factor 's'.
     */
    explicit DividedDifferences(size_t max_len = 10000, size_t s_max = 500)
        : max_len_(max_len),
          s_max_(s_max),
          N_(max_len + kExtraLength) {
        // Pre-allocate memory to avoid reallocations during simulation hot loops
        z_.reserve(max_len_);
        h_.resize(N_ + 1);
        g_matrix_.resize(max_len_ * s_max_);
        results_.resize(max_len_ + 1);
    }

    // The Rule of Zero: default copy/move/destructor are sufficient
    ~DividedDifferences() = default;
    DividedDifferences(const DividedDifferences&) = default;
    DividedDifferences& operator=(const DividedDifferences&) = default;
    DividedDifferences(DividedDifferences&&) = default;
    DividedDifferences& operator=(DividedDifferences&&) = default;

    /**
     * @brief Adds a new element to the input list and updates the DDEFs.
     * @param z_new The new input value to add.
     */
    void add_element(double z_new) {
        // Check if a full re-computation is needed due to buffer size or 's' change
        if (check_and_recompute_if_needed(z_new)) {
            return;
        }

        const size_t n = z_.size(); // Current size before adding
        z_.push_back(z_new);

        if (n == 0) { // First element initialization
            s_ = 1;
            mu_ = z_[0];
            exp_mu_ = ExExFloat::from_exp(mu_);

            h_[0] = ExExFloat(1.0);
            for (size_t k = 1; k <= N_; ++k) {
                h_[k] = h_[k - 1] / static_cast<double>(s_);
            }

            if (mu_ != z_[0]) {
                for (size_t k = N_; k > 0; --k) {
                    h_[k - 1] += h_[k] * ((z_[0] - mu_) / k);
                }
            }
            
            ExExFloat current_g = exp_mu_ * h_[0];
            for (size_t k = 0; k < s_ - 1; ++k) {
                g_matrix_[k * max_len_] = current_g;
                current_g *= h_[0];
            }
            results_[0] = current_g;

        } else { // Incremental update for subsequent elements
            for (size_t k = N_; k > n; --k) {
                h_[k - 1] += h_[k] * ((z_[n] - mu_) / k);
            }
            
            // Cache the value of h_[n] before the main transformation loop
            const ExExFloat h_n_val = h_[n];

            for (size_t k = n; k >= 1; --k) {
                h_[k - 1] = (h_[k - 1] * n + h_[k] * (z_[n] - z_[n - k])) / (n - k + 1.0);
            }

            ExExFloat current_g = exp_mu_ * h_n_val; // Use the cached value
            for (size_t k = 0; k < s_ - 1; ++k) {
                g_matrix_[k * max_len_ + n] = current_g;
                
                ExExFloat next_g = g_matrix_[k * max_len_] * h_[n];
                for (size_t j = 1; j <= n; ++j) {
                    next_g += g_matrix_[k * max_len_ + j] * h_[n - j];
                }
                current_g = next_g;
            }
            results_[n] = current_g;
        }
    }

    /**
     * @brief Removes the last element from the input list and updates DDEFs.
     */
    void remove_element() noexcept {
        if (z_.empty()) {
            return;
        }

        const size_t n = z_.size() - 1; // Size after removal
        const double z_removed = z_.back();
        z_.pop_back();

        if (empty()) { // If it was the last element, just clear state
            clear_state();
            return;
        }

        // Inverse transformation to restore the state of h_
        for (size_t k = 1; k <= n; ++k) {
            h_[k - 1] = (h_[k - 1] * (n - k + 1.0) - h_[k] * (z_removed - z_[n - k])) / n;
        }
        for (size_t k = n + 1; k <= N_; ++k) {
            h_[k - 1] -= h_[k] * ((z_removed - mu_) / k);
        }
    }

    /**
     * @brief Resets the calculator to its initial empty state.
     */
    void clear() noexcept {
        z_.clear();
        clear_state();
    }

    /**
     * @brief Gets the calculated divided difference exp[z_0, ..., z_index].
     * @param index The upper index of the divided difference.
     * @return The calculated DDEF value as an ExExFloat.
     */
    const ExExFloat& get_result(size_t index) const {
        if (index >= z_.size()) {
            throw std::out_of_range("Index out of range in get_result.");
        }
        return results_[index];
    }

    /// @brief Returns the number of input points.
    size_t size() const noexcept { return z_.size(); }
    /// @brief Checks if there are no input points.
    bool empty() const noexcept { return z_.empty(); }
    /// @brief Returns the current scaling factor 's'.
    int scaling_factor_s() const noexcept { return s_; }

private:
    // --- Constants from the paper ---
    static constexpr double kScalingDivisor = 3.5;
    static constexpr size_t kExtraLength = 30;

    // --- Internal State ---
    std::vector<double> z_;
    std::vector<ExExFloat> h_;
    std::vector<ExExFloat> g_matrix_; // 2D matrix (s_max x max_len) flattened
    std::vector<ExExFloat> results_;

    int s_ = 1;
    double mu_ = 0.0;
    ExExFloat exp_mu_;

    // --- Capacity ---
    size_t max_len_;
    size_t s_max_;
    size_t N_; // max_len_ + kExtraLength

    /**
     * @brief Resets internal calculation state without clearing inputs.
     */
    void clear_state() noexcept {
        s_ = 1;
        mu_ = 0.0;
        // No need to clear vectors, they will be overwritten.
    }

    /**
     * @brief Checks if a full re-computation is necessary and performs it.
     * @param z_new The next element to be added.
     * @return True if a re-computation was performed, false otherwise.
     */
    bool check_and_recompute_if_needed(double z_new) {
        if (z_.empty()) return false;

        const bool s_changed = std::abs(z_new - mu_) / kScalingDivisor > s_;
        const bool buffer_full = z_.size() >= max_len_;

        if (s_changed || buffer_full) {
            // Temporarily add the new element to calculate new s and mu
            z_.push_back(z_new);
            recompute_all();
            return true;
        }
        return false;
    }

    /**
     * @brief Performs a full re-computation of DDEFs from the current z_ vector.
     *
     * This is called when the scaling factor 's' changes or when internal
     * buffers need to be expanded.
     */
    void recompute_all() {
        if (z_.empty()) {
            clear_state();
            return;
        }

        // 1. Calculate new s and mu based on the full current z_ vector
        mu_ = std::accumulate(z_.begin(), z_.end(), 0.0) / z_.size();
        auto [min_z, max_z] = std::minmax_element(z_.begin(), z_.end());
        const double max_abs_diff = *max_z - *min_z;
        const int new_s = static_cast<int>(std::ceil(max_abs_diff / kScalingDivisor));
        
        // 2. Check if buffers need to be resized
        if (z_.size() >= max_len_ || new_s > s_max_) {
            size_t new_max_len = std::max(max_len_ * 2, z_.size());
            size_t new_s_max = std::max(s_max_ * 2, static_cast<size_t>(new_s));
            resize_internal_buffers(new_max_len, new_s_max);
        }

        // 3. Re-run the 'add_element' logic for all elements from scratch
        std::vector<double> z_copy = z_;
        clear(); // Fully reset the state
        
        for (const auto& val : z_copy) {
            add_element(val);
        }
    }

    /**
     * @brief Resizes all internal buffers to new capacities.
     */
    void resize_internal_buffers(size_t new_max_len, size_t new_s_max) {
        max_len_ = new_max_len;
        s_max_ = new_s_max;
        N_ = max_len_ + kExtraLength;

        z_.reserve(max_len_);
        h_.resize(N_ + 1);
        g_matrix_.resize(max_len_ * s_max_);
        results_.resize(max_len_ + 1);
    }
};

} // namespace pmrqmc::core