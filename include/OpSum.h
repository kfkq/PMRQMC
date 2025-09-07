#pragma once
#include <vector>
#include <complex>
#include <string>
#include <map>
#include <bitset>
#include <variant>
#include <optional>

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

} // namespace pmrqmc