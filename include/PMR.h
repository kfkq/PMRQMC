#pragma once
#include <vector>
#include <complex>
#include <string>
#include <map>
#include <bitset>
#include <variant>
#include <optional>

namespace pmrqmc {

// Forward declaration for OpSum
class OpSum;

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
    void add(const OpSum& opsum);
    
    // Access individual OpSum objects
    OpSum& operator[](size_t index);
    const OpSum& operator[](size_t index) const;
    
    // Get all OpSum objects
    const std::vector<OpSum>& get_opsums() const;
    
    // Get the number of OpSum objects in the bulk
    size_t size() const;
    
    // Check if bulk is empty
    bool empty() const;
    
    // Clear all OpSum objects
    void clear();
    
    // Get the name of the bulk collection
    const std::string& get_name() const;
    
    // Set the name of the bulk collection
    void set_name(const std::string& name);
    
    // Iterator support for range-based for loops
    std::vector<OpSum>::iterator begin();
    std::vector<OpSum>::iterator end();
    std::vector<OpSum>::const_iterator begin() const;
    std::vector<OpSum>::const_iterator end() const;
    std::vector<OpSum>::const_iterator cbegin() const;
    std::vector<OpSum>::const_iterator cend() const;
};

// Observable Builder class
class ObservableBuilder {
public:
    // Create single observable from OpSum and Hamiltonian PMR
    static PMRResult create_single(const OpSum& observable, const PMRResult& hamiltonian_pmr);
    
    // Create bulk observable from OpSumBulk and Hamiltonian PMR
    static PMRResult create_bulk(const OpSumBulk& observables, const PMRResult& hamiltonian_pmr);
};

// Main PMR function
PMRResult pmr(const OpSum& hamiltonian);

// Observable PMR functions
PMRResult pmr_observable(const OpSum& observable, const PMRResult& hamiltonian_pmr);
PMRResult pmr_observable_bulk(const OpSumBulk& observables, const PMRResult& hamiltonian_pmr);

// Helper functions
std::vector<bool> int_to_bitset(int value, int size);
std::vector<bool> pauli_string_to_bitset(const std::vector<std::pair<std::string, int>>& pauli_ops, int n_qubits);
std::complex<double> pauli_coefficient(const std::vector<std::pair<std::string, int>>& pauli_ops);

} // namespace pmrqmc