#pragma once
#include <vector>
#include <complex>
#include <string>
#include <map>
#include <bitset>
#include <variant>

namespace pmrqmc {

class OpSum {
private:
    struct OperatorTerm {
        std::complex<double> coefficient;
        std::vector<std::pair<std::string, int>> pauli_ops; // ('X', 0), ('Y', 1), ('Sz', 2), etc.
    };
    
    std::vector<OperatorTerm> terms;
    int max_qubit_index = -1;
    std::string name_;
    
public:
    OpSum() = default;
    explicit OpSum(const std::string& name) : name_(name) {}
    
    // Set a name for the Hamiltonian (for debugging/identification)
    void set_name(const std::string& name) { name_ = name; }
    const std::string& get_name() const { return name_; }
    
    // Add operator term with vector of Pauli operators
    void add(double coefficient, const std::vector<std::pair<std::string, int>>& pauli_ops);
    void add(std::complex<double> coefficient, const std::vector<std::pair<std::string, int>>& pauli_ops);
    
    // Add operator term with alternating Pauli chars and qubit indices
    void add(double coefficient, std::initializer_list<std::pair<std::string, int>> pauli_ops);
    void add(std::complex<double> coefficient, std::initializer_list<std::pair<std::string, int>> pauli_ops);
    
    // Flexible add method - handles multiple formats (Pauli chars and qubit ints)
    template<typename... Args>
    void add(double coefficient, const std::string& first_op, Args&&... args) {
        std::vector<std::pair<std::string, int>> pauli_ops;
        add_term_impl(pauli_ops, first_op, std::forward<Args>(args)...);
        add(coefficient, pauli_ops);
    }
    
    template<typename... Args>
    void add(std::complex<double> coefficient, const std::string& first_op, Args&&... args) {
        std::vector<std::pair<std::string, int>> pauli_ops;
        add_term_impl(pauli_ops, first_op, std::forward<Args>(args)...);
        add(coefficient, pauli_ops);
    }
    
    // Getters
    const std::vector<OperatorTerm>& get_terms() const { return terms; }
    int get_max_qubit_index() const { return max_qubit_index + 1; } // Number of qubits
    
    // Clear all terms
    void clear() { 
        terms.clear(); 
        max_qubit_index = -1; 
    }
    
    // Print Hamiltonian (for debugging)
    void print() const;
    
private:
    // Helper for template argument parsing
    void add_term_impl(std::vector<std::pair<std::string, int>>& pauli_ops);
    
    // Helper functions
    bool is_valid_pauli_operator(const std::string& op);
    std::vector<std::pair<std::string, int>> convert_spin_operators(const std::vector<std::pair<std::string, int>>& ops);
    
    template<typename First, typename... Rest>
    void add_term_impl(std::vector<std::pair<std::string, int>>& pauli_ops, 
                      First&& first, Rest&&... rest) {
        if constexpr (std::is_same_v<std::decay_t<First>, std::string> || 
                      std::is_same_v<std::decay_t<First>, const char*> ||
                      std::is_same_v<std::decay_t<First>, char>) {
            // Pauli operator - push it with a placeholder
            std::string pauli_str;
            if constexpr (std::is_same_v<std::decay_t<First>, char>) {
                pauli_str = std::string(1, first);
            } else {
                pauli_str = first;
            }
            pauli_ops.emplace_back(pauli_str, -1); // -1 as placeholder
            add_term_impl(pauli_ops, std::forward<Rest>(rest)...);
        } else if constexpr (std::is_same_v<std::decay_t<First>, int>) {
            // Qubit index - assign to the last Pauli operator
            int qubit = first;
            if (!pauli_ops.empty() && pauli_ops.back().second == -1) {
                pauli_ops.back().second = qubit;
            } else {
                throw std::invalid_argument("Qubit index must follow Pauli operator");
            }
            add_term_impl(pauli_ops, std::forward<Rest>(rest)...);
        } else {
            throw std::invalid_argument("Arguments must be Pauli operator (string/char) followed by int (qubit)");
        }
    }
};

// PMR result structure - directly usable by QMC code
struct PMRResult {
    int N;                          // Number of qubits
    int Nop;                        // Number of permutation operators
    int Ncycles;                    // Number of fundamental cycles
    
    // Permutation matrices and cycles - directly usable
    std::vector<std::vector<bool>> P_matrix;
    std::vector<std::vector<bool>> cycles;
    
    // Diagonal terms (D0)
    std::vector<std::complex<double>> D0_coeff;
    std::vector<std::vector<bool>> D0_product;
    int D0_size = 0;
    
    // Off-diagonal terms (D)
    std::vector<std::vector<std::complex<double>>> D_coeff;
    std::vector<std::vector<std::vector<bool>>> D_product;
    std::vector<int> D_size;
    int D_maxsize = 0;
    
    // Helper methods for direct use in QMC
    bool has_diagonal_terms() const { return D0_size > 0; }
    bool has_offdiagonal_terms() const { return Nop > 0; }
    
    // Get permutation matrix for operator i
    const std::vector<bool>& get_permutation(int i) const { 
        return P_matrix[i]; 
    }
    
    // Get diagonal coefficient for j-th diagonal term
    std::complex<double> get_diagonal_coeff(int j) const { 
        return D0_coeff[j]; 
    }
    
    // Get off-diagonal coefficients for operator i
    const std::vector<std::complex<double>>& get_offdiagonal_coeffs(int i) const { 
        return D_coeff[i]; 
    }
    
    // Get number of terms for operator i
    int get_operator_terms(int i) const { 
        return D_size[i]; 
    }
};

// Main PMR function
PMRResult pmr(const OpSum& hamiltonian);

// Helper functions
std::vector<bool> int_to_bitset(int value, int size);
std::vector<bool> pauli_string_to_bitset(const std::vector<std::pair<std::string, int>>& pauli_ops, int n_qubits);
std::complex<double> pauli_coefficient(const std::vector<std::pair<std::string, int>>& pauli_ops);

} // namespace pmrqmc