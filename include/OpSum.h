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

// Observable PMR result structure
struct ObservablePMRResult {
    // Observable mapping to Hamiltonian permutation operators
    std::vector<int> MNop;           // Number of permutation operators for each observable
    std::vector<std::vector<std::vector<bool>>> MP;  // Mapping matrix for each observable
    
    // Diagonal terms for observables (MD0)
    std::vector<int> MD0_size;       // Number of diagonal terms for each observable
    std::vector<std::vector<std::complex<double>>> MD0_coeff;  // Diagonal coefficients
    std::vector<std::vector<std::vector<bool>>> MD0_product;   // Diagonal products
    
    // Off-diagonal terms for observables (MD)
    std::vector<std::vector<int>> MD_size;       // Number of off-diagonal terms for each observable
    std::vector<std::vector<std::vector<std::complex<double>>>> MD_coeff;  // Off-diagonal coefficients
    std::vector<std::vector<std::vector<std::vector<bool>>>> MD_product;  // Off-diagonal products
    std::vector<int> MD_maxsize;     // Maximum number of off-diagonal terms for each observable
    
    // Observable names
    std::vector<std::string> Mnames;
    
    // Helper methods
    int get_num_observables() const { return MNop.size(); }
    
    // Check if observable has diagonal terms
    bool has_diagonal_terms(int obs_idx) const { 
        return obs_idx < MD0_size.size() && MD0_size[obs_idx] > 0; 
    }
    
    // Check if observable has off-diagonal terms
    bool has_offdiagonal_terms(int obs_idx) const { 
        return obs_idx < MNop.size() && MNop[obs_idx] > 0; 
    }
    
    // Get diagonal coefficients for observable
    const std::vector<std::complex<double>>& get_diagonal_coeffs(int obs_idx) const { 
        return MD0_coeff[obs_idx]; 
    }
    
    // Get off-diagonal coefficients for observable
    const std::vector<std::vector<std::complex<double>>>& get_offdiagonal_coeffs(int obs_idx) const { 
        return MD_coeff[obs_idx]; 
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
    
    // Observable PMR data (if computed for observables)
    std::optional<ObservablePMRResult> observable_data;
    
    // Helper methods for direct use in QMC
    bool has_diagonal_terms() const { return D0_size > 0; }
    bool has_offdiagonal_terms() const { return Nop > 0; }
    bool has_observable_data() const { return observable_data.has_value(); }
    
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
    
    // Get observable data
    const ObservablePMRResult& get_observable_data() const { 
        return observable_data.value(); 
    }
};

// Generic collection of OpSum objects for bulk observables
class OpSumBulk {
private:
    std::vector<OpSum> opsums;
    std::string name_;
    
public:
    explicit OpSumBulk(const std::string& name = "bulk") : name_(name) {}
    
    // Add an OpSum to the bulk collection
    void add(const OpSum& opsum) {
        opsums.push_back(opsum);
    }
    
    // Access individual OpSum objects
    OpSum& operator[](size_t index) {
        return opsums[index];
    }
    
    const OpSum& operator[](size_t index) const {
        return opsums[index];
    }
    
    // Get all OpSum objects
    const std::vector<OpSum>& get_opsums() const {
        return opsums;
    }
    
    // Get the number of OpSum objects in the bulk
    size_t size() const {
        return opsums.size();
    }
    
    // Check if bulk is empty
    bool empty() const {
        return opsums.empty();
    }
    
    // Clear all OpSum objects
    void clear() {
        opsums.clear();
    }
    
    // Get the name of the bulk collection
    const std::string& get_name() const {
        return name_;
    }
    
    // Set the name of the bulk collection
    void set_name(const std::string& name) {
        name_ = name;
    }
    
    // Iterator support for range-based for loops
    auto begin() { return opsums.begin(); }
    auto end() { return opsums.end(); }
    auto begin() const { return opsums.begin(); }
    auto end() const { return opsums.end(); }
    auto cbegin() const { return opsums.cbegin(); }
    auto cend() const { return opsums.cend(); }
};

// Main PMR function
PMRResult pmr(const OpSum& hamiltonian);

// Observable PMR functions
PMRResult pmr_observable(const OpSum& observable, const PMRResult& hamiltonian_pmr);
PMRResult pmr_observable_bulk(const OpSumBulk& observables, const PMRResult& hamiltonian_pmr);

// Observable Builder class
class ObservableBuilder {
public:
    // Create single observable from OpSum and Hamiltonian PMR
    static PMRResult create_single(const OpSum& observable, const PMRResult& hamiltonian_pmr) {
        return pmr_observable(observable, hamiltonian_pmr);
    }
    
    // Create bulk observable from OpSumBulk and Hamiltonian PMR
    static PMRResult create_bulk(const OpSumBulk& observables, const PMRResult& hamiltonian_pmr) {
        return pmr_observable_bulk(observables, hamiltonian_pmr);
    }
};

// Helper functions
std::vector<bool> int_to_bitset(int value, int size);
std::vector<bool> pauli_string_to_bitset(const std::vector<std::pair<std::string, int>>& pauli_ops, int n_qubits);
std::complex<double> pauli_coefficient(const std::vector<std::pair<std::string, int>>& pauli_ops);

} // namespace pmrqmc