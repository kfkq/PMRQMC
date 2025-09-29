#include <divdiff.hpp>
#include <stdexcept>
#include <cassert>
#include <array>

namespace pmrqmc {

namespace {

// Global precomputed inverse powers of 2 for performance
std::unique_ptr<std::array<double, ExExFloat::max_exponent>> inv_powers_of_two;

} // anonymous namespace

// Global initialization function
void divdiff_init() {
    if (!inv_powers_of_two) {
        inv_powers_of_two = std::make_unique<std::array<double, ExExFloat::max_exponent>>();
        double current = 1.0;
        for (size_t i = 0; i < ExExFloat::max_exponent; ++i) {
            (*inv_powers_of_two)[i] = current;
            current /= 2.0;
        }
    }
}

void divdiff_cleanup() noexcept {
    inv_powers_of_two.reset();
}

// ==================== DivDiff::Builder Implementation ====================

DivDiff::Builder& DivDiff::Builder::max_length(int n) {
    if (n <= 0) {
        throw std::invalid_argument("max_length must be positive");
    }
    max_len_ = n;
    return *this;
}

DivDiff::Builder& DivDiff::Builder::max_scaling_factor(int s) {
    if (s <= 0) {
        throw std::invalid_argument("max_scaling_factor must be positive");
    }
    max_scaling_ = s;
    return *this;
}

DivDiff DivDiff::Builder::build() const {
    return DivDiff(max_len_, max_scaling_);
}

std::unique_ptr<DivDiff> DivDiff::Builder::build_unique() const noexcept {
    try {
        return std::make_unique<DivDiff>(max_len_, max_scaling_);
    } catch (...) {
        return nullptr;
    }
}

// ==================== DivDiff Implementation ====================

DivDiff::DivDiff(int max_length, int max_scaling_factor)
    : max_length_(max_length)
    , max_scaling_(max_scaling_factor) {
    if (max_length <= 0 || max_scaling_factor <= 0) {
        throw std::invalid_argument("DivDiff dimensions must be positive");
    }
    allocate_memory();
}

void DivDiff::add_element(double energy) {
    if (current_length_ >= max_length_) {
        // Reinitialize with all current elements
        add_all_elements(current_length_);
        return;
    }

    energies_.push_back(energy);

    if (current_length_ == 0) {
        initialize_single_point();
    } else if (needs_scaling_change()) {
        add_all_elements(current_length_ + 1);
    } else {
        compute_new_element();
    }

    ++current_length_;
}

void DivDiff::remove_element() {
    if (current_length_ <= 0) {
        return;
    }

    // Implementation based on legacy algorithm
    int n = current_length_ - 1;
    int N = max_length_ + 10;  // extra length

    for (int k = 1; k <= n; ++k) {
        h_matrix_[k - 1] = (h_matrix_[k - 1] * n + h_matrix_[k] * (energies_[n] - energies_[n - k])) / (n - k + 1);
    }

    for (int k = n + 1; k <= N; ++k) {
        h_matrix_[k - 1] -= h_matrix_[k] * (energies_[n] - mean_energy_);
    }

    energies_.pop_back();
    --current_length_;
}

bool DivDiff::remove_value(double value) {
    auto it = std::find(energies_.begin(), energies_.end(), value);
    if (it == energies_.end()) {
        return false;
    }

    int index = std::distance(energies_.begin(), it);
    int n = current_length_ - 1;

    // Remove elements before the target
    for (int j = n; j >= index; --j) {
        remove_element();
    }

    // Add back elements after the target
    for (int j = index; j < n; ++j) {
        add_element(energies_[j + 1]);
    }

    return true;
}

const std::vector<ExExFloat>& DivDiff::divdiffs() const noexcept {
    return divdiffs_;
}

const std::vector<double>& DivDiff::energies() const noexcept {
    return energies_;
}

int DivDiff::current_length() const noexcept {
    return current_length_;
}

int DivDiff::scaling_factor() const noexcept {
    return scaling_factor_;
}

size_t DivDiff::memory_usage() const noexcept {
    return sizeof(DivDiff) +
           energies_.capacity() * sizeof(double) +
           h_matrix_.capacity() * sizeof(ExExFloat) +
           divdiffs_.capacity() * sizeof(ExExFloat) +
           ddd_matrix_.capacity() * sizeof(ExExFloat);
}

void DivDiff::print_state(std::ostream& os) const {
    os << "DivDiff State:\n";
    os << "  Current length: " << current_length_ << "\n";
    os << "  Scaling factor: " << scaling_factor_ << "\n";
    os << "  Mean energy: " << mean_energy_ << "\n";
    os << "  Energies: [";
    for (size_t i = 0; i < energies_.size(); ++i) {
        os << energies_[i];
        if (i < energies_.size() - 1) os << ",";
    }
    os << "]\n";
    os << "  DivDiffs: [";
    for (size_t i = 0; i < divdiffs_.size(); ++i) {
        divdiffs_[i].print(os);
        if (i < divdiffs_.size() - 1) os << ",";
    }
    os << "]\n";
}

void DivDiff::initialize_single_point() {
    scaling_factor_ = 1;
    mean_energy_ = energies_[0];

    const int N = max_length_ + 10;
    const int maxlen = max_length_;

    // Initialize h matrix: h[0] = 1, h[k] = h[k-1]/s for k=1 to N
    h_matrix_[0] = ExExFloat(1.0);
    for (int k = 1; k <= N; ++k) {
        h_matrix_[k] = h_matrix_[k - 1] / static_cast<double>(scaling_factor_);
    }

    // Adjust for shifted mean energy (if mu != z[0])
    if (mean_energy_ != energies_[0]) {
        for (int k = N; k > 0; --k) {
            ExExFloat mu_diff = energies_[0] - mean_energy_;
            h_matrix_[k - 1] += h_matrix_[k] * (mu_diff / static_cast<double>(k));
        }
    }

    // Compute exp(z[0]) and initialize ddd matrix
    ExExFloat exp_z0(0.0);
    exp_z0.init_expmu(energies_[0]);  // exp(z[0])

    ExExFloat curr = exp_z0 * h_matrix_[0];
    for (int k = 0; k < scaling_factor_ - 1; ++k) {
        ddd_matrix_[k * maxlen] = curr;
        curr *= h_matrix_[0];
    }

    // divdiffs[0] = exp(z[0]) (the DDEF f[z0] = exp(z0))
    divdiffs_[0].init_expmu(energies_[0]);
}

void DivDiff::add_all_elements(int len) {
    current_length_ = 0;

    if (len <= 0) return;

    // Calculate mean and max diff for scaling
    double sum = 0.0;
    for (int i = 0; i < len; ++i) {
        sum += energies_[i];
    }
    double mu = sum / len;

    double max_diff{};
    for (int i = 0; i < len; ++i) {
        max_diff = std::max(max_diff, std::abs(energies_[i] - mu));
    }

    scaling_factor_ = std::max(1, static_cast<int>(std::ceil(max_diff / 3.5)));
    mean_energy_ = mu;

    // Resize if needed
    resize_if_needed();

    // Add each element with new scaling
    for (int i = 0; i < len; ++i) {
        add_element(energies_[i]);
    }
}

bool DivDiff::needs_scaling_change() const noexcept {
    if (current_length_ == 0) return false;
    double latest_energy = energies_.back();
    return std::abs(latest_energy - mean_energy_) / 3.5 > scaling_factor_;
}

void DivDiff::compute_new_element() {
    // Full implementation of the DDEF (Divided Differences of Exponential Function) algorithm
    // Based on the Gupta et al. paper and legacy AddElement implementation

    int n = current_length_;
    int N = static_cast<int>(h_matrix_.size()) - 1;  // Total matrix size
    int maxlen = max_length_;

    // Create exponentiated mean energy
    ExExFloat expmu(mean_energy_);

    // Update h matrix: for(k=N;k>n;k--) h[k-1] += h[k]*(z[n]-mu)/k
    for (int k = N; k > n; --k) {
        ExExFloat mu_diff = energies_[n] - mean_energy_;
        ExExFloat update = h_matrix_[k] * (mu_diff / static_cast<double>(k));
        h_matrix_[k - 1] += update;
    }

    // Compute expmu*h[n]
    ExExFloat curr = expmu * h_matrix_[n];

    // Update h matrix for the new element: for(k=n;k>=1;k--) h[k-1] = (h[k-1]*n + h[k]*(z[n]-z[n-k]))/(n-k+1);
    for (int k = n; k >= 1; --k) {
        ExExFloat mu_diff2 = energies_[n] - energies_[n - k];
        ExExFloat update2 = h_matrix_[k] * mu_diff2;
        h_matrix_[k - 1] = (h_matrix_[k - 1] * static_cast<double>(n) + update2) / static_cast<double>(n - k + 1);
    }

    // Update ddd matrix and compute current ddd value:
    // for(k=0;k<s-1;k++){
    //     ddd[k*maxlen+n] = curr;
    //     curr = ddd[k*maxlen]*h[n]; for(j=1;j<=n;j++) curr += ddd[k*maxlen+j]*h[n-j];
    // }
    for (int k = 0; k < scaling_factor_ - 1; ++k) {
        // Store current curr in ddd matrix
        ddd_matrix_[k * maxlen + n] = curr;

        // Compute next curr: ddd[k*maxlen]*h[n] + sum_{j=1 to n} ddd[k*maxlen+j]*h[n-j]
        curr = ddd_matrix_[k * maxlen] * h_matrix_[n];

        for (int j = 1; j <= n; ++j) {
            curr += ddd_matrix_[k * maxlen + j] * h_matrix_[n - j];
        }
    }

    // Final result: divdiffs[n] = curr
    divdiffs_[n] = curr;
}


void DivDiff::allocate_memory() {
    int total_size = max_length_ + 10 + 1;  // extra length + 1 for safety

    try {
        energies_.reserve(max_length_);
        h_matrix_.resize(total_size);
        divdiffs_.resize(total_size);
        ddd_matrix_.resize(max_length_ * max_scaling_);
    } catch (const std::bad_alloc&) {
        free_memory();
        throw std::runtime_error("Failed to allocate memory for DivDiff");
    }
}

void DivDiff::free_memory() noexcept {
    energies_.clear();
    h_matrix_.clear();
    divdiffs_.clear();
    ddd_matrix_.clear();
}

void DivDiff::resize_if_needed() {
    // Check if current allocations are sufficient
    if (energies_.capacity() < static_cast<size_t>(max_length_)) {
        energies_.reserve(max_length_);
    }

    int total_size = max_length_ + 10 + 1;
    if (h_matrix_.size() < static_cast<size_t>(total_size)) {
        h_matrix_.resize(total_size);
        divdiffs_.resize(total_size);
    }

    if (ddd_matrix_.size() < static_cast<size_t>(max_length_ * max_scaling_)) {
        ddd_matrix_.resize(max_length_ * max_scaling_);
    }
}

// ==================== ExExFloat Implementation ====================

ExExFloat::ExExFloat(double val) noexcept : mantissa(val), exponent(0) {
    normalize();
}

ExExFloat& ExExFloat::operator=(double rhs) noexcept {
    mantissa = rhs;
    exponent = 0;
    normalize();
    return *this;
}

void ExExFloat::normalize() noexcept {
    if (mantissa == 0.0) {
        exponent = 0;
        return;
    }

    int temp_exponent;
    mantissa = std::frexp(mantissa, &temp_exponent);
    exponent += temp_exponent;
}

ExExFloat ExExFloat::operator+(const ExExFloat& rhs) const noexcept {
    ExExFloat result;

    if (rhs.exponent >= exponent) {
        if (rhs.mantissa != 0.0) {
            result.mantissa = rhs.mantissa + mantissa * (*inv_powers_of_two)[rhs.exponent - exponent];
            result.exponent = rhs.exponent;
            result.normalize();
        } else {
            result.mantissa = mantissa;
            result.exponent = exponent;
        }
    } else {
        if (mantissa != 0.0) {
            result.mantissa = mantissa + rhs.mantissa * (*inv_powers_of_two)[exponent - rhs.exponent];
            result.exponent = exponent;
            result.normalize();
        } else {
            result.mantissa = rhs.mantissa;
            result.exponent = rhs.exponent;
        }
    }

    return result;
}

ExExFloat ExExFloat::operator-(const ExExFloat& rhs) const noexcept {
    ExExFloat result;

    if (rhs.exponent >= exponent) {
        if (rhs.mantissa != 0.0) {
            result.mantissa = mantissa * (*inv_powers_of_two)[rhs.exponent - exponent] - rhs.mantissa;
            result.exponent = rhs.exponent;
            result.normalize();
        } else {
            result.mantissa = mantissa;
            result.exponent = exponent;
        }
    } else {
        if (mantissa != 0.0) {
            result.mantissa = mantissa - rhs.mantissa * (*inv_powers_of_two)[exponent - rhs.exponent];
            result.exponent = exponent;
            result.normalize();
        } else {
            result.mantissa = -rhs.mantissa;
            result.exponent = rhs.exponent;
        }
    }

    return result;
}

ExExFloat ExExFloat::operator*(const ExExFloat& rhs) const noexcept {
    ExExFloat result;
    result.mantissa = mantissa * rhs.mantissa;
    result.exponent = exponent + rhs.exponent;
    result.normalize();
    return result;
}

ExExFloat ExExFloat::operator/(const ExExFloat& rhs) const {
    if (rhs.mantissa == 0.0) {
        throw std::domain_error("Division by zero in ExExFloat");
    }
    ExExFloat result;
    result.mantissa = mantissa / rhs.mantissa;
    result.exponent = exponent - rhs.exponent;
    result.normalize();
    return result;
}

ExExFloat ExExFloat::operator*(double rhs) const noexcept {
    ExExFloat result;
    result.mantissa = mantissa * rhs;
    result.exponent = exponent;
    result.normalize();
    return result;
}

ExExFloat ExExFloat::operator/(double rhs) const {
    if (rhs == 0.0) {
        throw std::domain_error("Division by zero in ExExFloat");
    }
    ExExFloat result;
    result.mantissa = mantissa / rhs;
    result.exponent = exponent;
    result.normalize();
    return result;
}

ExExFloat operator*(double lhs, const ExExFloat& rhs) noexcept {
    return rhs * lhs;
}

ExExFloat operator/(double lhs, const ExExFloat& rhs) {
    ExExFloat result;
    if (rhs.mantissa == 0.0) {
        throw std::domain_error("Division by zero in ExExFloat");
    }
    result.mantissa = lhs / rhs.mantissa;
    result.exponent = -rhs.exponent;
    result.normalize();
    return result;
}

ExExFloat& ExExFloat::operator+=(const ExExFloat& rhs) noexcept {
    *this = *this + rhs;
    return *this;
}

ExExFloat& ExExFloat::operator-=(const ExExFloat& rhs) noexcept {
    *this = *this - rhs;
    return *this;
}

ExExFloat& ExExFloat::operator*=(const ExExFloat& rhs) noexcept {
    *this = *this * rhs;
    return *this;
}

ExExFloat& ExExFloat::operator/=(const ExExFloat& rhs) {
    *this = *this / rhs;
    return *this;
}

ExExFloat& ExExFloat::operator*=(double rhs) noexcept {
    *this = *this * rhs;
    return *this;
}

ExExFloat& ExExFloat::operator/=(double rhs) {
    *this = *this / rhs;
    return *this;
}

bool ExExFloat::operator>=(const ExExFloat& rhs) const noexcept {
    if (mantissa == 0.0) return rhs.mantissa <= 0.0;
    if (rhs.mantissa == 0.0) return mantissa >= 0.0;

    if (exponent > rhs.exponent) return true;
    if (exponent < rhs.exponent) return false;
    return mantissa >= rhs.mantissa;
}

bool ExExFloat::operator>=(double rhs) const noexcept {
    if (rhs == 0.0) return mantissa >= 0.0;

    ExExFloat r(rhs);
    return *this >= r;
}

double ExExFloat::double_value() const noexcept {
    return std::ldexp(mantissa, exponent);
}

int ExExFloat::sign() const noexcept {
    return (mantissa > 0.0) - (mantissa < 0.0);
}

ExExFloat ExExFloat::abs() const noexcept {
    ExExFloat result;
    result.mantissa = std::fabs(mantissa);
    result.exponent = exponent;
    return result;
}

ExExFloat ExExFloat::sqrt() const {
    if (mantissa < 0.0) {
        throw std::domain_error("Square root of negative number in ExExFloat");
    }

    ExExFloat result;
    if (exponent % 2 == 0) {
        result.mantissa = std::sqrt(mantissa);
        result.exponent = exponent / 2;
    } else {
        result.mantissa = std::sqrt(2.0 * mantissa);
        result.exponent = (exponent - 1) / 2;
    }
    result.normalize();
    return result;
}

void ExExFloat::init_expmu(double mu) {
    double log10_e = 0.4342944819032518;  // log10(e)
    double exp_part = mu * log10_e;
    exponent = static_cast<int>(std::ceil(exp_part));
    mantissa = std::pow(2.0, exp_part - std::ceil(exp_part));
    normalize();
}

void ExExFloat::print(std::ostream& os) const {
    double double_val = double_value();
    if (exponent < 50 && exponent > -50) {
        os << double_val;
    } else {
        // For very large numbers, print in scientific notation
        double exp10 = exponent * 0.30102999566398114; // log10(2)
        double m = mantissa * std::pow(10.0, exp10 - std::floor(exp10));
        double exp10_floor = std::floor(exp10);
        if (std::fabs(m) < 1.0) {
            exp10_floor--;
            m *= 10.0;
        }
        os << m << "e" << exp10_floor;
    }
}

} // namespace pmrqmc