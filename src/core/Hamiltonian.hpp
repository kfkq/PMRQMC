#pragma once

#include "OperatorTerm.hpp"
#include "ExExFloat.hpp"
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace pmrqmc::core {

class Hamiltonian {
private:
    std::vector<OperatorTerm> terms_;
    std::vector<OperatorTerm> diagonal_terms_;
    std::vector<OperatorTerm> off_diagonal_terms_;
    int num_qubits_ = 0;
    bool sorted_ = false;

public:
    Hamiltonian() = default;
    
    // Constructor from vector of terms
    explicit Hamiltonian(const std::vector<OperatorTerm>& terms) : terms_(terms) {
        compute_num_qubits();
        separate_diagonal_off_diagonal();
        sort_terms();
    }
    
    // Factory function to load from file
    static Hamiltonian from_file(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        std::vector<OperatorTerm> terms;
        std::string line;
        
        while (std::getline(file, line)) {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#') {
                continue;
            }
            
            std::istringstream iss(line);
            std::complex<double> coeff;
            std::string pauli_str;
            
            // Parse coefficient (can be real or complex)
            std::string coeff_str;
            iss >> coeff_str;
            
            // Check if coefficient is complex
            if (coeff_str.find('i') != std::string::npos || 
                coeff_str.find('j') != std::string::npos) {
                // Complex coefficient parsing
                size_t plus_pos = coeff_str.find('+');
                size_t minus_pos = coeff_str.find('-', 1); // Skip first character
                
                double real_part = 0.0;
                double imag_part = 0.0;
                
                if (plus_pos != std::string::npos) {
                    // Format like "1.0+2.0i" or "1.0 + 2.0i"
                    std::string real_str = coeff_str.substr(0, plus_pos);
                    std::string imag_str = coeff_str.substr(plus_pos + 1);
                    
                    real_part = std::stod(real_str);
                    imag_part = std::stod(imag_str);
                } else if (minus_pos != std::string::npos) {
                    // Format like "1.0-2.0i" or "1.0 - 2.0i"
                    std::string real_str = coeff_str.substr(0, minus_pos);
                    std::string imag_str = coeff_str.substr(minus_pos);
                    
                    real_part = std::stod(real_str);
                    imag_part = std::stod(imag_str);
                } else {
                    // Format like "2.0i" or "-3.0i"
                    imag_part = std::stod(coeff_str);
                }
                
                coeff = std::complex<double>(real_part, imag_part);
            } else {
                // Real coefficient
                coeff = std::complex<double>(std::stod(coeff_str), 0.0);
            }
            
            // Parse Pauli string
            iss >> pauli_str;
            
            // Create and add the term
            terms.emplace_back(coeff, pauli_str);
        }
        
        file.close();
        
        return Hamiltonian(terms);
    }
    
    // Factory function to load from legacy format (like in the original code)
    static Hamiltonian from_legacy_format(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        std::vector<OperatorTerm> terms;
        std::string line;
        
        while (std::getline(file, line)) {
            // Skip empty lines
            if (line.empty()) continue;
            
            std::istringstream iss(line);
            std::string token;
            
            // Parse the line format: coefficient qubit1 op1 qubit2 op2 ...
            std::vector<std::string> tokens;
            while (iss >> token) {
                tokens.push_back(token);
            }
            
            if (tokens.size() < 2 || tokens.size() % 2 != 1) {
                continue; // Invalid format, skip
            }
            
            // Parse coefficient
            double coeff_real = 0.0;
            double coeff_imag = 0.0;
            
            // Check if first token is complex
            if (tokens[0].find('i') != std::string::npos || 
                tokens[0].find('j') != std::string::npos) {
                // Complex coefficient
                coeff_imag = std::stod(tokens[0].substr(0, tokens[0].length() - 1));
            } else {
                // Real coefficient
                coeff_real = std::stod(tokens[0]);
            }
            
            // Build Pauli string
            std::string pauli_str;
            for (size_t i = 1; i < tokens.size(); i += 2) {
                int qubit = std::stoi(tokens[i]);
                std::string op_str = tokens[i + 1];
                
                // Normalize to uppercase
                std::transform(op_str.begin(), op_str.end(), op_str.begin(), ::toupper);
                
                if (op_str != "I") {  // Only add non-identity operators
                    pauli_str += op_str + std::to_string(qubit);
                }
            }
            
            if (pauli_str.empty()) {
                pauli_str = "I";  // Identity operator
            }
            
            terms.emplace_back(std::complex<double>(coeff_real, coeff_imag), pauli_str);
        }
        
        file.close();
        return Hamiltonian(terms);
    }
    
    // Accessors
    int num_qubits() const { return num_qubits_; }
    size_t num_terms() const { return terms_.size(); }
    size_t num_diagonal_terms() const { return diagonal_terms_.size(); }
    size_t num_off_diagonal_terms() const { return off_diagonal_terms_.size(); }
    
    const std::vector<OperatorTerm>& get_terms() const { return terms_; }
    const std::vector<OperatorTerm>& get_diagonal_terms() const { return diagonal_terms_; }
    const std::vector<OperatorTerm>& get_off_diagonal_terms() const { return off_diagonal_terms_; }
    
    std::vector<OperatorTerm>& get_terms() { return terms_; }
    std::vector<OperatorTerm>& get_diagonal_terms() { return diagonal_terms_; }
    std::vector<OperatorTerm>& get_off_diagonal_terms() { return off_diagonal_terms_; }
    
    // Add a term to the Hamiltonian
    void add_term(const OperatorTerm& term) {
        terms_.push_back(term);
        if (term.is_diagonal()) {
            diagonal_terms_.push_back(term);
        } else {
            off_diagonal_terms_.push_back(term);
        }
        update_num_qubits(term);
        sorted_ = false;
    }
    
    // Remove terms with very small coefficients (numerical cleanup)
    void cleanup(double threshold = 1e-15) {
        auto remove_predicate = [threshold](const OperatorTerm& term) {
            return std::abs(term.coefficient()) < threshold;
        };
        
        terms_.erase(std::remove_if(terms_.begin(), terms_.end(), remove_predicate), terms_.end());
        diagonal_terms_.erase(std::remove_if(diagonal_terms_.begin(), diagonal_terms_.end(), remove_predicate), diagonal_terms_.end());
        off_diagonal_terms_.erase(std::remove_if(off_diagonal_terms_.begin(), off_diagonal_terms_.end(), remove_predicate), off_diagonal_terms_.end());
    }
    
    // Check if Hamiltonian is Hermitian (with numerical tolerance)
    bool is_hermitian(double tolerance = 1e-15) const {
        for (const auto& term : terms_) {
            bool found_hermitian_pair = false;
            
            if (std::abs(term.coefficient().imag()) < tolerance) {
                // Real coefficient - term should be its own Hermitian conjugate
                found_hermitian_pair = true;
            } else {
                // Complex coefficient - look for Hermitian conjugate
                std::complex<double> herm_coeff = std::conj(term.coefficient());
                for (const auto& other_term : terms_) {
                    if (term.same_pauli_string(other_term) && 
                        std::abs(other_term.coefficient() - herm_coeff) < tolerance) {
                        found_hermitian_pair = true;
                        break;
                    }
                }
            }
            
            if (!found_hermitian_pair) {
                return false;
            }
        }
        
        return true;
    }
    
    // Get maximum coefficient magnitude (for numerical analysis)
    double max_coefficient_magnitude() const {
        double max_mag = 0.0;
        for (const auto& term : terms_) {
            double mag = std::abs(term.coefficient());
            max_mag = std::max(max_mag, mag);
        }
        return max_mag;
    }
    
    // Get minimum coefficient magnitude (excluding zeros)
    double min_coefficient_magnitude() const {
        double min_mag = std::numeric_limits<double>::max();
        bool found_nonzero = false;
        for (const auto& term : terms_) {
            double mag = std::abs(term.coefficient());
            if (mag > 0.0) {  // Only exclude actual zeros
                min_mag = std::min(min_mag, mag);
                found_nonzero = true;
            }
        }
        return found_nonzero ? min_mag : 0.0;
    }
    
    // Sort terms by coefficient magnitude (useful for analysis)
    void sort_terms() {
        std::sort(terms_.begin(), terms_.end(), 
            [](const OperatorTerm& a, const OperatorTerm& b) {
                return std::abs(a.coefficient()) > std::abs(b.coefficient());
            });
        
        std::sort(diagonal_terms_.begin(), diagonal_terms_.end(), 
            [](const OperatorTerm& a, const OperatorTerm& b) {
                return std::abs(a.coefficient()) > std::abs(b.coefficient());
            });
        
        std::sort(off_diagonal_terms_.begin(), off_diagonal_terms_.end(), 
            [](const OperatorTerm& a, const OperatorTerm& b) {
                return std::abs(a.coefficient()) > std::abs(b.coefficient());
            });
        
        sorted_ = true;
    }
    
    // Convert to legacy format for compatibility with original code
    std::string to_legacy_format() const {
        std::ostringstream oss;
        
        for (const auto& term : terms_) {
            auto coeff = term.coefficient();
            oss << coeff.real();
            if (std::abs(coeff.imag()) > 1e-15) {
                oss << (coeff.imag() >= 0 ? "+" : "") << coeff.imag() << "i";
            }
            
            const auto& pauli = term.pauli_string();
            for (int qubit : pauli.get_qubits()) {
                PauliOp op = pauli.get_operator(qubit);
                oss << " " << qubit << " " << to_string(op);
            }
            
            oss << "\n";
        }
        
        return oss.str();
    }
    
    // Calculate energy expectation value for a given state
    template<typename StateType>
    std::complex<double> energy_expectation(const StateType& state) const {
        std::complex<double> energy = 0.0;
        
        for (const auto& term : terms_) {
            auto [new_state, phase, coeff] = term.apply(state);
            // This is a simplified version - real QMC would need proper state handling
            energy += coeff * phase;
        }
        
        return energy;
    }
    
    // Check if Hamiltonian has only real coefficients
    bool is_real() const {
        for (const auto& term : terms_) {
            if (std::abs(term.coefficient().imag()) > 1e-15) {
                return false;
            }
        }
        return true;
    }
    
    // Check if Hamiltonian has no off-diagonal terms (classical Hamiltonian)
    bool is_diagonal() const {
        return off_diagonal_terms_.empty();
    }
    
    // Print summary
    void print_summary(std::ostream& os = std::cout) const {
        os << "Hamiltonian Summary:\n";
        os << "  Number of qubits: " << num_qubits_ << "\n";
        os << "  Total terms: " << terms_.size() << "\n";
        os << "  Diagonal terms: " << diagonal_terms_.size() << "\n";
        os << "  Off-diagonal terms: " << off_diagonal_terms_.size() << "\n";
        os << "  Hermitian: " << (is_hermitian() ? "Yes" : "No") << "\n";
        os << "  Real coefficients: " << (is_real() ? "Yes" : "No") << "\n";
        os << "  Max coefficient magnitude: " << max_coefficient_magnitude() << "\n";
        os << "  Min coefficient magnitude: " << min_coefficient_magnitude() << "\n";
    }

private:
    // Compute number of qubits from the terms
    void compute_num_qubits() {
        num_qubits_ = 0;
        for (const auto& term : terms_) {
            update_num_qubits(term);
        }
    }
    
    // Update number of qubits based on a term
    void update_num_qubits(const OperatorTerm& term) {
        int term_max_qubit = term.max_qubit();
        num_qubits_ = std::max(num_qubits_, term_max_qubit + 1);  // +1 for 0-based indexing
    }
    
    // Separate terms into diagonal and off-diagonal
    void separate_diagonal_off_diagonal() {
        diagonal_terms_.clear();
        off_diagonal_terms_.clear();
        
        for (const auto& term : terms_) {
            if (term.is_diagonal()) {
                diagonal_terms_.push_back(term);
            } else {
                off_diagonal_terms_.push_back(term);
            }
        }
    }
};

}