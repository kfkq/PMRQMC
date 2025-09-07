#include <OpSum.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace pmrqmc {

// OpSum constructor removed - string parameter was unnecessary

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
    
    // Add the term with converted operators and final coefficient
    terms.push_back({final_coeff, converted_ops});
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
    
    // Add the term with converted operators and final coefficient
    terms.push_back({final_coeff, converted_ops});
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

// PMR computation functions
std::vector<bool> int_to_bitset(int value, int size) {
    std::vector<bool> result(size, false);
    for (int i = 0; i < size && value > 0; ++i) {
        result[i] = (value & (1 << i));
    }
    return result;
}

std::vector<bool> pauli_string_to_bitset(const std::vector<std::pair<std::string, int>>& pauli_ops, int n_qubits) {
    std::vector<bool> result(n_qubits, false);
    
    for (const auto& [pauli, qubit] : pauli_ops) {
        if (pauli == "X" || pauli == "Y" || pauli == "Sx" || pauli == "Sy" || pauli == "S+" || pauli == "S-") {
            if (qubit < n_qubits) {
                result[qubit] = true;
            }
        }
        // 'Z', 'Sz', 'I' don't contribute to permutation bitset
    }
    
    return result;
}

std::complex<double> pauli_coefficient(const std::vector<std::pair<std::string, int>>& pauli_ops) {
    std::complex<double> coeff(1.0, 0.0);
    
    for (const auto& [pauli, qubit] : pauli_ops) {
        if (pauli == "Y") {
            coeff *= std::complex<double>(0.0, 1.0); // Y = i * (|0><1| - |1><0|)
        }
    }
    
    return coeff;
}

// GF(2) vector addition
std::vector<int> gf2_add(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors must be same size for GF(2) addition");
    }
    
    std::vector<int> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = (vec1[i] + vec2[i]) % 2;
    }
    return result;
}

// Compute nullspace (adapted from legacy code)
std::vector<std::vector<int>> compute_nullspace(const std::vector<std::vector<int>>& matrix) {
    if (matrix.empty()) return {};
    
    int numRows = matrix.size();
    int numCols = matrix[0].size();
    std::vector<std::vector<int>> matrix_RE = matrix;
    std::vector<int> marked_rows;
    
    // Gaussian elimination in GF(2)
    for (int j = 0; j < numCols; ++j) {
        for (int i = 0; i < numRows; ++i) {
            if (matrix_RE[i][j] == 1) {
                marked_rows.push_back(i);
                for (int k = 0; k < numCols; ++k) {
                    if (k != j && matrix_RE[i][k] == 1) {
                        for (int l = 0; l < numRows; ++l) {
                            matrix_RE[l][k] = matrix_RE[l][j] ^ matrix_RE[l][k];
                        }
                    }
                }
                break;
            }
        }
    }
    
    // Remove duplicates and sort
    std::sort(marked_rows.begin(), marked_rows.end());
    marked_rows.erase(std::unique(marked_rows.begin(), marked_rows.end()), marked_rows.end());
    
    std::vector<std::vector<int>> marked_matrix;
    for (int row : marked_rows) {
        marked_matrix.push_back(matrix_RE[row]);
    }
    
    std::vector<std::vector<int>> nullspace;
    for (int i = 0; i < numRows; ++i) {
        auto it = std::find(marked_rows.begin(), marked_rows.end(), i);
        if (it == marked_rows.end()) {
            std::vector<int> nullvector(numRows, 0);
            nullvector[i] = 1;
            
            std::vector<int> row_i = matrix_RE[i];
            for (int j = 0; j < numCols; ++j) {
                if (row_i[j] == 1) {
                    for (int k = 0; k < marked_matrix.size(); ++k) {
                        if (marked_matrix[k][j] == 1) {
                            nullvector[marked_rows[k]] = 1;
                        }
                    }
                }
            }
            nullspace.push_back(nullvector);
        }
    }
    
    return nullspace;
}

// Minimize cycle sizes (adapted from legacy code)
int minimize_cycles(std::vector<std::vector<int>>& null_eigs) {
    int changes_made = 0;
    int nullsize = null_eigs.size();
    
    // Sort by number of 1s
    std::sort(null_eigs.begin(), null_eigs.end(), 
              [](const std::vector<int>& a, const std::vector<int>& b) {
                  int sum_a = 0, sum_b = 0;
                  for (int x : a) sum_a += x;
                  for (int x : b) sum_b += x;
                  return sum_a < sum_b;
              });
    
    for (int k = nullsize - 1; k > 0; --k) {
        int null_k = 0;
        for (int x : null_eigs[k]) null_k += x;
        
        for (int m = 0; m < k; ++m) {
            std::vector<int> curr = gf2_add(null_eigs[k], null_eigs[m]);
            int curr_sum = 0;
            for (int x : curr) curr_sum += x;
            
            if (curr_sum < null_k) {
                null_eigs[k] = curr;
                changes_made = 1;
                break;
            }
        }
    }
    
    return changes_made;
}

PMRResult pmr(const OpSum& hamiltonian) {
    PMRResult result;
    
    // Get basic information
    const auto& terms = hamiltonian.get_terms();
    result.N = hamiltonian.get_max_qubit_index();
    
    if (terms.empty()) {
        return result; // Empty Hamiltonian
    }
    
    // Group terms by their permutation representation
    std::map<std::vector<bool>, std::vector<std::pair<std::complex<double>, std::vector<bool>>>> permutation_groups;
    
    for (const auto& term : terms) {
        auto permutation = pauli_string_to_bitset(term.pauli_ops, result.N);
        auto coeff = term.coefficient * pauli_coefficient(term.pauli_ops);
        
        // Convert Z operators to bitset representation
        std::vector<bool> z_product(result.N, false);
        for (const auto& [pauli, qubit] : term.pauli_ops) {
            if (pauli == "Z" || pauli == "Y" || pauli == "Sz") {
                if (qubit < result.N) {
                    z_product[qubit] = !z_product[qubit];
                }
            }
        }
        
        permutation_groups[permutation].emplace_back(coeff, z_product);
    }
    
    // Remove zero coefficients
    for (auto it = permutation_groups.begin(); it != permutation_groups.end(); ) {
        auto& group = it->second;
        group.erase(
            std::remove_if(group.begin(), group.end(), 
                          [](const auto& pair) { return std::abs(pair.first) < 1e-10; }),
            group.end()
        );
        if (group.empty()) {
            it = permutation_groups.erase(it);
        } else {
            ++it;
        }
    }
    
    // Sort by permutation (like legacy code)
    std::vector<std::pair<std::vector<bool>, decltype(permutation_groups)::mapped_type>> sorted_groups;
    for (auto& [perm, group] : permutation_groups) {
        sorted_groups.emplace_back(perm, std::move(group));
    }
    
    std::sort(sorted_groups.begin(), sorted_groups.end(),
              [](const auto& a, const auto& b) {
                  return a.first < b.first;
              });
    
    // Separate diagonal (identity permutation) and off-diagonal terms
    bool has_diagonal = false;
    std::vector<std::vector<bool>> nontrivial_permutations;
    std::vector<decltype(sorted_groups)::value_type> offdiagonal_groups;
    
    for (const auto& [perm, group] : sorted_groups) {
        bool is_identity = std::all_of(perm.begin(), perm.end(), [](bool x) { return !x; });
        
        if (is_identity) {
            has_diagonal = true;
            // Handle diagonal terms (D0)
            for (const auto& [coeff, z_prod] : group) {
                result.D0_coeff.push_back(coeff);
                result.D0_product.push_back(z_prod);
            }
            result.D0_size = result.D0_coeff.size();
        } else {
            nontrivial_permutations.push_back(perm);
            offdiagonal_groups.emplace_back(perm, group);
        }
    }
    
    result.Nop = nontrivial_permutations.size();
    
    // Convert permutation matrices to integer format for nullspace computation
    std::vector<std::vector<int>> perm_int;
    for (const auto& perm : nontrivial_permutations) {
        std::vector<int> perm_row;
        for (bool bit : perm) {
            perm_row.push_back(bit ? 1 : 0);
        }
        perm_int.push_back(perm_row);
    }
    
    // Compute cycles (nullspace)
    if (!perm_int.empty()) {
        auto nullspace = compute_nullspace(perm_int);
        while (minimize_cycles(nullspace));
        
        result.Ncycles = nullspace.size();
        for (const auto& cycle : nullspace) {
            result.cycles.push_back(std::vector<bool>(result.Nop, false));
            for (int i = 0; i < result.Nop; ++i) {
                result.cycles.back()[i] = (cycle[i] == 1);
            }
        }
    } else {
        result.Ncycles = 0;
    }
    
    // Store permutation matrices
    result.P_matrix = std::move(nontrivial_permutations);
    
    // Handle off-diagonal terms
    result.D_maxsize = 0;
    for (const auto& [perm, group] : offdiagonal_groups) {
        std::vector<std::complex<double>> coeffs;
        std::vector<std::vector<bool>> products;
        
        for (const auto& [coeff, z_prod] : group) {
            coeffs.push_back(coeff);
            products.push_back(z_prod);
        }
        
        result.D_coeff.push_back(coeffs);
        result.D_product.push_back(products);
        result.D_size.push_back(coeffs.size());
        result.D_maxsize = std::max(result.D_maxsize, (int)coeffs.size());
    }
    
    // Ensure consistent sizing for D_coeff
    for (auto& coeffs : result.D_coeff) {
        coeffs.resize(result.D_maxsize, 0.0);
    }
    
    return result;
}

} // namespace pmrqmc