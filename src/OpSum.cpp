#include <OpSum.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace pmrqmc {

void OpSum::add(double coefficient, const std::vector<std::pair<std::string, int>>& pauli_ops) {
    if (pauli_ops.empty()) {
        throw std::invalid_argument("At least one Pauli operator required");
    }
    
    // Validate Pauli operators and find max qubit index
    for (const auto& [pauli, qubit] : pauli_ops) {
        if (!is_valid_pauli_operator(pauli)) {
            throw std::invalid_argument("Pauli operator must be 'X', 'Y', 'Z', 'I', 'Sx', 'Sy', 'Sz', 'S+', or 'S-'");
        }
        if (qubit < 0) {
            throw std::invalid_argument("Qubit index must be non-negative");
        }
        max_qubit_index = std::max(max_qubit_index, qubit);
    }
    
    // Handle spin operators by converting them to Pauli operators with proper coefficients
    auto converted_ops = convert_spin_operators(pauli_ops);
    
    // Apply spin operator coefficient scaling
    double final_coeff = coefficient;
    for (const auto& [op, qubit] : pauli_ops) {
        if (op == "Sx" || op == "Sy" || op == "Sz") {
            final_coeff *= 0.5;  // S = σ/2
        }
    }

    // Apply Pauli coefficient for phase (e.g., Y contributing i)
    std::complex<double> phase_coeff(1.0, 0.0);
    for (const auto& [pauli, qubit] : converted_ops) {
        if (pauli == "Y") {
            phase_coeff *= std::complex<double>(0.0, 1.0);
        }
    }
    // Add the term with converted operators and final coefficient
    terms.push_back({final_coeff * phase_coeff, converted_ops});
}

void OpSum::add(double coefficient, std::initializer_list<std::pair<std::string, int>> pauli_ops) {
    std::vector<std::pair<std::string, int>> ops_vec(pauli_ops.begin(), pauli_ops.end());
    add(coefficient, ops_vec);
}

void OpSum::add(std::complex<double> coefficient, const std::vector<std::pair<std::string, int>>& pauli_ops) {
    if (pauli_ops.empty()) {
        throw std::invalid_argument("At least one Pauli operator required");
    }
    
    // Validate Pauli operators and find max qubit index
    for (const auto& [pauli, qubit] : pauli_ops) {
        if (!is_valid_pauli_operator(pauli)) {
            throw std::invalid_argument("Pauli operator must be 'X', 'Y', 'Z', 'I', 'Sx', 'Sy', 'Sz', 'S+', or 'S-'");
        }
        if (qubit < 0) {
            throw std::invalid_argument("Qubit index must be non-negative");
        }
        max_qubit_index = std::max(max_qubit_index, qubit);
    }
    
    // Handle spin operators by converting them to Pauli operators with proper coefficients
    auto converted_ops = convert_spin_operators(pauli_ops);
    
    // Apply spin operator coefficient scaling
    std::complex<double> final_coeff = coefficient;
    for (const auto& [op, qubit] : pauli_ops) {
        if (op == "Sx" || op == "Sy" || op == "Sz") {
            final_coeff *= 0.5;  // S = σ/2
        }
    }

    // Apply Pauli coefficient for phase (e.g., Y contributing i)
    std::complex<double> phase_coeff(1.0, 0.0);
    for (const auto& [pauli, qubit] : converted_ops) {
        if (pauli == "Y") {
            phase_coeff *= std::complex<double>(0.0, 1.0);
        }
    }
    // Add the term with converted operators and final coefficient
    terms.push_back({final_coeff * phase_coeff, converted_ops});
}

void OpSum::add(std::complex<double> coefficient, std::initializer_list<std::pair<std::string, int>> pauli_ops) {
    std::vector<std::pair<std::string, int>> ops_vec(pauli_ops.begin(), pauli_ops.end());
    add(coefficient, ops_vec);
}

// Helper function to validate Pauli operators
bool OpSum::is_valid_pauli_operator(const std::string& op) {
    return op == "X" || op == "Y" || op == "Z" || op == "I" ||
           op == "Sx" || op == "Sy" || op == "Sz";
}

// Helper function to convert spin operators to Pauli operators
std::vector<std::pair<std::string, int>> OpSum::convert_spin_operators(const std::vector<std::pair<std::string, int>>& ops) {
    std::vector<std::pair<std::string, int>> result;
    
    for (const auto& [op, qubit] : ops) {
        if (op == "Sx") {
            result.emplace_back("X", qubit);
        } else if (op == "Sy") {
            result.emplace_back("Y", qubit);
        } else if (op == "Sz") {
            result.emplace_back("Z", qubit);
        } else {
            result.emplace_back(op, qubit);
        }
    }
    
    return result;
}

// Base case for template recursion
void OpSum::add_term_impl(std::vector<std::pair<std::string, int>>& pauli_ops) {
    // Base case - no more arguments
}

void OpSum::print() const {
    std::cout << "Hamiltonian with " << terms.size() << " terms:" << std::endl;
    for (const auto& term : terms) {
        std::cout << term.coefficient << " * ";
        for (const auto& [pauli, qubit] : term.pauli_ops) {
            std::cout << pauli << "_" << qubit << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace pmrqmc