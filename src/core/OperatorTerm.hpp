#pragma once

#include "Pauli.hpp"
#include <complex>

namespace pmrqmc::core {

// Represents a single term in a Hamiltonian: coefficient * PauliString
class OperatorTerm {
private:
    std::complex<double> coefficient_;
    PauliString pauli_string_;

public:
    OperatorTerm() : coefficient_(0.0) {}
    
    OperatorTerm(std::complex<double> coeff, const PauliString& pauli)
        : coefficient_(coeff), pauli_string_(pauli) {}
    
    OperatorTerm(std::complex<double> coeff, const std::string& pauli_str)
        : coefficient_(coeff), pauli_string_(PauliString(pauli_str)) {}
    
    // Accessors
    const std::complex<double>& coefficient() const { return coefficient_; }
    std::complex<double>& coefficient() { return coefficient_; }
    
    const PauliString& pauli_string() const { return pauli_string_; }
    PauliString& pauli_string() { return pauli_string_; }
    
    // Check if this term is diagonal (only Z and I operators)
    bool is_diagonal() const {
        for (int qubit : pauli_string_.get_qubits()) {
            PauliOp op = pauli_string_.get_operator(qubit);
            if (op == PauliOp::X || op == PauliOp::Y) {
                return false;
            }
        }
        return true;
    }
    
    // Check if this term is off-diagonal (contains X or Y operators)
    bool is_off_diagonal() const {
        return !is_diagonal();
    }
    
    // Get the maximum qubit index involved
    int max_qubit() const {
        std::vector<int> qubits = pauli_string_.get_qubits();
        if (qubits.empty()) return 0;
        return *std::max_element(qubits.begin(), qubits.end());
    }
    
    // Check if this term acts on a specific qubit
    bool acts_on_qubit(int qubit) const {
        return pauli_string_.has_operator(qubit);
    }
    
    // Count the number of qubits this term acts on
    int num_qubits() const {
        return static_cast<int>(pauli_string_.size());
    }
    
    // Check if this term is the identity
    bool is_identity() const {
        return pauli_string_.is_identity();
    }
    
    // Scalar multiplication
    OperatorTerm operator*(std::complex<double> scalar) const {
        return OperatorTerm(coefficient_ * scalar, pauli_string_);
    }
    
    // In-place scalar multiplication
    OperatorTerm& operator*=(std::complex<double> scalar) {
        coefficient_ *= scalar;
        return *this;
    }
    
    // Term equality (ignores coefficient differences)
    bool same_pauli_string(const OperatorTerm& other) const {
        return pauli_string_ == other.pauli_string_;
    }
    
    // Full equality
    bool operator==(const OperatorTerm& other) const {
        return coefficient_ == other.coefficient_ && pauli_string_ == other.pauli_string_;
    }
    
    bool operator!=(const OperatorTerm& other) const {
        return !(*this == other);
    }
    
    // Convert to string representation
    std::string to_string() const {
        std::string result;
        
        // Format coefficient
        if (coefficient_.imag() == 0.0) {
            if (coefficient_.real() == 1.0 && !pauli_string_.is_identity()) {
                // Omit coefficient 1 for non-identity terms
            } else if (coefficient_.real() == -1.0 && !pauli_string_.is_identity()) {
                result += "-";
            } else {
                result += std::to_string(coefficient_.real()) + " * ";
            }
        } else {
            result += "(" + std::to_string(coefficient_.real()) + " + " + 
                     std::to_string(coefficient_.imag()) + "i) * ";
        }
        
        result += pauli_string_.to_string();
        return result;
    }
    
    // Apply to a basis state
    // Returns the resulting state, phase factor, and coefficient contribution
    std::tuple<std::string, std::complex<double>, std::complex<double>> 
    apply(const std::string& state) const {
        auto [new_state, phase] = pauli_string_.apply(state);
        return {new_state, phase, coefficient_};
    }
};

// Free functions for OperatorTerm arithmetic
inline OperatorTerm operator*(std::complex<double> scalar, const OperatorTerm& term) {
    return term * scalar;
}

inline OperatorTerm operator*(double scalar, const OperatorTerm& term) {
    return term * std::complex<double>(scalar, 0.0);
}

inline OperatorTerm operator*(const OperatorTerm& term, double scalar) {
    return term * std::complex<double>(scalar, 0.0);
}

}