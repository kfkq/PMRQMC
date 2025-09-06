#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <compare>
#include <array>

// This class is a modernized C++ implementation based on the concepts from:
// L. Gupta, L. Barash, I. Hen, Computer Physics Communications 254, 107385 (2020).
//
// It provides an extended-exponent floating-point number type to avoid overflow/underflow
// in standard double-precision arithmetic during QMC calculations.

namespace pmrqmc::core {

class ExExFloat {
private:
    double mantissa_ = 0.0;
    int exponent_ = 0; // Represents the number: mantissa_ * 2^exponent_

    // Cache for powers of 2 - initialized once per thread
    struct PowerCache {
        std::array<double, 128> powers;
        bool initialized = false;
        
        PowerCache() {
            // Initialize cache on first use
            for (int i = 0; i < 128; ++i) {
                powers[i] = std::pow(2.0, -i);
            }
            initialized = true;
        }
    };
    
    static thread_local PowerCache cache_;

    void normalize() noexcept {
        if (mantissa_ == 0.0 || std::isinf(mantissa_) || std::isnan(mantissa_)) {
            exponent_ = 0;
            return;
        }
        int exp;
        mantissa_ = std::frexp(mantissa_, &exp);
        exponent_ += exp;
    }

public:
    // --- Constructors ---
    constexpr ExExFloat() noexcept = default;

    explicit ExExFloat(double val) noexcept : mantissa_(val) {
        normalize();
    }

    // --- Factory Functions ---
    // Creates an ExExFloat representing exp(mu)
    static ExExFloat from_exp(double mu) noexcept {
        ExExFloat res;
        // Using the identity: exp(mu) = 2^(mu * log2(e))
        constexpr double log2_e = 1.4426950408889634;
        double total_exp = mu * log2_e;
        res.exponent_ = static_cast<int>(std::ceil(total_exp));
        res.mantissa_ = std::pow(2.0, total_exp - res.exponent_);
        // The mantissa will be in (0.5, 1.0], so no need to normalize.
        return res;
    }

    // --- Conversion and Accessors ---
    double to_double() const noexcept {
        return std::ldexp(mantissa_, exponent_);
    }

    constexpr int sgn() const noexcept {
        return (mantissa_ > 0.0) - (mantissa_ < 0.0);
    }

    // --- Unary Operators ---
    ExExFloat operator-() const noexcept {
        ExExFloat res = *this;
        res.mantissa_ = -res.mantissa_;
        return res;
    }

    // --- Mathematical Functions ---
    ExExFloat abs() const noexcept {
        ExExFloat res = *this;
        res.mantissa_ = std::fabs(res.mantissa_);
        return res;
    }

    ExExFloat sqrt() const noexcept {
        if (mantissa_ < 0.0) {
            return ExExFloat(std::numeric_limits<double>::quiet_NaN());
        }
        if (mantissa_ == 0.0) {
            return ExExFloat(0.0);
        }

        ExExFloat res;
        if (exponent_ % 2 == 0) {
            res.mantissa_ = std::sqrt(mantissa_);
            res.exponent_ = exponent_ / 2;
        } else {
            res.mantissa_ = std::sqrt(2.0 * mantissa_);
            res.exponent_ = (exponent_ - 1) / 2;
        }
        res.normalize();
        return res;
    }

    // Additional mathematical functions for QMC applications
    ExExFloat square() const noexcept {
        ExExFloat res = *this;
        res.mantissa_ *= res.mantissa_;
        res.exponent_ *= 2;
        res.normalize();
        return res;
    }

    bool is_zero() const noexcept {
        return mantissa_ == 0.0;
    }

    bool is_finite() const noexcept {
        return std::isfinite(mantissa_);
    }

    bool is_nan() const noexcept {
        return std::isnan(mantissa_);
    }

    bool is_inf() const noexcept {
        return std::isinf(mantissa_);
    }

    // --- Compound Assignment Operators ---
    ExExFloat& operator+=(const ExExFloat& rhs) noexcept {
        if (rhs.mantissa_ == 0.0) {
            return *this;
        }
        if (mantissa_ == 0.0) {
            *this = rhs;
            return *this;
        }

        if (exponent_ >= rhs.exponent_) {
            int exp_diff = exponent_ - rhs.exponent_;
            if (exp_diff < 128) {
                mantissa_ += rhs.mantissa_ * cache_.powers[exp_diff];
            } else {
                mantissa_ += std::ldexp(rhs.mantissa_, rhs.exponent_ - exponent_);
            }
        } else {
            int exp_diff = rhs.exponent_ - exponent_;
            if (exp_diff < 128) {
                mantissa_ = rhs.mantissa_ + mantissa_ * cache_.powers[exp_diff];
                exponent_ = rhs.exponent_;
            } else {
                mantissa_ = rhs.mantissa_ + std::ldexp(mantissa_, exponent_ - rhs.exponent_);
                exponent_ = rhs.exponent_;
            }
        }
        normalize();
        return *this;
    }

    ExExFloat& operator-=(const ExExFloat& rhs) noexcept {
        if (rhs.mantissa_ == 0.0) {
            return *this;
        }
        if (mantissa_ == 0.0) {
            *this = -rhs;
            return *this;
        }

        if (exponent_ >= rhs.exponent_) {
            int exp_diff = exponent_ - rhs.exponent_;
            if (exp_diff < 128) {
                mantissa_ -= rhs.mantissa_ * cache_.powers[exp_diff];
            } else {
                mantissa_ -= std::ldexp(rhs.mantissa_, rhs.exponent_ - exponent_);
            }
        } else {
            int exp_diff = rhs.exponent_ - exponent_;
            if (exp_diff < 128) {
                mantissa_ = std::ldexp(mantissa_, exponent_ - rhs.exponent_) - rhs.mantissa_;
                exponent_ = rhs.exponent_;
            } else {
                mantissa_ = std::ldexp(mantissa_, exponent_ - rhs.exponent_) - rhs.mantissa_;
                exponent_ = rhs.exponent_;
            }
        }
        normalize();
        return *this;
    }

    ExExFloat& operator*=(const ExExFloat& rhs) noexcept {
        mantissa_ *= rhs.mantissa_;
        exponent_ += rhs.exponent_;
        normalize();
        return *this;
    }

    ExExFloat& operator*=(double rhs) noexcept {
        mantissa_ *= rhs;
        normalize();
        return *this;
    }

    ExExFloat& operator/=(const ExExFloat& rhs) noexcept {
        mantissa_ /= rhs.mantissa_;
        exponent_ -= rhs.exponent_;
        normalize();
        return *this;
    }

    ExExFloat& operator/=(double rhs) noexcept {
        mantissa_ /= rhs;
        normalize();
        return *this;
    }

    // --- Comparison Operator (C++20 Spaceship Operator) ---
    std::partial_ordering operator<=>(const ExExFloat& rhs) const noexcept {
        if (mantissa_ == 0.0 && rhs.mantissa_ == 0.0) {
            return std::partial_ordering::equivalent;
        }
        int s1 = sgn();
        int s2 = rhs.sgn();
        if (s1 != s2) {
            return s1 <=> s2;
        }
        // Signs are the same from here on
        if (s1 == 0) { // Both are zero
            return std::partial_ordering::equivalent;
        }
        if (s1 > 0) { // Both are positive
            if (exponent_ != rhs.exponent_) {
                return exponent_ <=> rhs.exponent_;
            }
            return mantissa_ <=> rhs.mantissa_;
        }
        // Both are negative
        if (exponent_ != rhs.exponent_) {
            return rhs.exponent_ <=> exponent_; // Reversed
        }
        return rhs.mantissa_ <=> mantissa_; // Reversed
    }
    
    // We need an explicit operator== for cases where <=> is not sufficient (e.g., unordered maps)
    bool operator==(const ExExFloat& rhs) const noexcept {
        return (exponent_ == rhs.exponent_ && mantissa_ == rhs.mantissa_) ||
               (mantissa_ == 0.0 && rhs.mantissa_ == 0.0);
    }
};

// Define the static thread-local cache member
thread_local ExExFloat::PowerCache ExExFloat::cache_;

// --- Free Functions for Binary and Mixed-Type Arithmetic ---

inline ExExFloat operator+(ExExFloat lhs, const ExExFloat& rhs) noexcept {
    lhs += rhs;
    return lhs;
}
inline ExExFloat operator-(ExExFloat lhs, const ExExFloat& rhs) noexcept {
    lhs -= rhs;
    return lhs;
}
inline ExExFloat operator*(ExExFloat lhs, const ExExFloat& rhs) noexcept {
    lhs *= rhs;
    return lhs;
}
inline ExExFloat operator/(ExExFloat lhs, const ExExFloat& rhs) noexcept {
    lhs /= rhs;
    return lhs;
}

// Mixed-type with double
inline ExExFloat operator+(ExExFloat lhs, double rhs) noexcept { return lhs + ExExFloat(rhs); }
inline ExExFloat operator+(double lhs, const ExExFloat& rhs) noexcept { return ExExFloat(lhs) + rhs; }
inline ExExFloat operator-(ExExFloat lhs, double rhs) noexcept { return lhs - ExExFloat(rhs); }
inline ExExFloat operator-(double lhs, const ExExFloat& rhs) noexcept { return ExExFloat(lhs) - rhs; }
inline ExExFloat operator*(ExExFloat lhs, double rhs) noexcept { return lhs * ExExFloat(rhs); }
inline ExExFloat operator*(double lhs, const ExExFloat& rhs) noexcept { return ExExFloat(lhs) * rhs; }
inline ExExFloat operator/(ExExFloat lhs, double rhs) noexcept { return lhs / ExExFloat(rhs); }
inline ExExFloat operator/(double lhs, const ExExFloat& rhs) noexcept { return ExExFloat(lhs) / rhs; }


// --- Stream Output ---

inline std::ostream& operator<<(std::ostream& os, const ExExFloat& val) {
    if (val.sgn() == 0) {
        os << "0.0";
        return os;
    }
    // For human-readable output, convert to base 10 representation
    // val = m * 2^e = m_10 * 10^e_10
    // log10(val) = log10(m) + e * log10(2)
    constexpr double log10_2 = 0.3010299956639812;
    double exp10 = std::log10(std::fabs(val.to_double()));
    double mant10 = val.to_double() / std::pow(10.0, std::floor(exp10));
    
    os << mant10 << "e" << std::floor(exp10);
    return os;
}

} // namespace pmrqmc::core