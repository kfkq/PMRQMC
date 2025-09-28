#include <PMR.h>
#include <OpSum.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace pmrqmc {

// Helper functions
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
                result[n_qubits - 1 - qubit] = true;  // REVERSED: to match legacy bit ordering
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

// Compute nullspace (adapted from legacy code) - FIXED VERSION
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
            std::vector<int> nullvector(numRows, 0);  // Fixed: nullvector should have numRows elements
            nullvector[i] = 1;
            
            std::vector<int> row_i = matrix_RE[i];
            for (int j = 0; j < numCols; ++j) {
                if (row_i[j] == 1) {
                    for (size_t k = 0; k < marked_matrix.size(); ++k) {
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

// Main PMR computation functions
PMR pmr(const OpSum& hamiltonian) {
    PMR result;
    
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
                    z_product[result.N - 1 - qubit] = !z_product[result.N - 1 - qubit];  // REVERSED: to match legacy bit ordering
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

// Observable PMR computation functions
PMR pmr_obs(const OpSum& observable, const PMR& hamiltonian_pmr) {
    OpSumBulk obs_bulk("single");
    obs_bulk.add(observable);
    return pmr_obs(obs_bulk, hamiltonian_pmr);
}

PMR pmr_obs(const OpSumBulk& observables, const PMR& hamiltonian_pmr) {
    PMR result = hamiltonian_pmr; // Copy Hamiltonian PMR data
    
    // Initialize observable data structure
    ObservablePMR obs_data;
    
    // Process each observable in the bulk
    for (const auto& observable : observables.get_opsums()) {
        obs_data.Mnames.push_back(observable.get_name());
        
        // Get observable terms
        const auto& terms = observable.get_terms();
        
        // Initialize observable-specific arrays
        obs_data.MNop.push_back(0);
        obs_data.MD0_size.push_back(0);
        obs_data.MD_maxsize.push_back(0);
        
        // Group observable terms by Hamiltonian permutation operators
        std::map<std::vector<bool>, std::vector<std::pair<std::complex<double>, std::vector<bool>>>> obs_groups;
        
        for (const auto& term : terms) {
            auto permutation = pauli_string_to_bitset(term.pauli_ops, result.N);
            auto coeff = term.coefficient * pauli_coefficient(term.pauli_ops);
            
            // Convert Z operators to bitset representation
            std::vector<bool> z_product(result.N, false);
            for (const auto& [pauli, qubit] : term.pauli_ops) {
                if (pauli == "Z" || pauli == "Y" || pauli == "Sz") {
                    if (qubit < result.N) {
                        z_product[result.N - 1 - qubit] = !z_product[result.N - 1 - qubit];  // REVERSED: to match legacy bit ordering
                    }
                }
            }
            
            obs_groups[permutation].emplace_back(coeff, z_product);
        }
        
        // Remove zero coefficients
        for (auto it = obs_groups.begin(); it != obs_groups.end(); ) {
            auto& group = it->second;
            group.erase(
                std::remove_if(group.begin(), group.end(), 
                              [](const auto& pair) { return std::abs(pair.first) < 1e-10; }),
                group.end()
            );
            if (group.empty()) {
                it = obs_groups.erase(it);
            } else {
                ++it;
            }
        }
        
        // Check if observable can be expressed in Hamiltonian permutation basis
        // Using nullspace method as in legacy prepare.cpp
        for (const auto& [obs_perm, group] : obs_groups) {
            // Skip identity (diagonal) operators - they're always valid
            bool is_identity = std::all_of(obs_perm.begin(), obs_perm.end(), [](bool x) { return !x; });
            if (is_identity) {
                continue;
            }
            
            // Convert Hamiltonian permutations to integer matrix for nullspace computation
            std::vector<std::vector<int>> ham_perm_int;
            for (int i = 0; i < result.Nop; ++i) {
                std::vector<int> perm_row;
                for (bool bit : result.P_matrix[i]) {
                    perm_row.push_back(bit ? 1 : 0);
                }
                ham_perm_int.push_back(perm_row);
            }
            
            // Convert observable permutation to integer vector
            std::vector<int> obs_perm_int;
            for (bool bit : obs_perm) {
                obs_perm_int.push_back(bit ? 1 : 0);
            }
            
            // Create combined matrix for nullspace computation
            std::vector<std::vector<int>> combined_matrix = ham_perm_int;
            combined_matrix.push_back(obs_perm_int);
            
            // Compute nullspace
            auto nullspace = compute_nullspace(combined_matrix);
            
            // Check if any nullspace vector has the last element set to 1
            // This means the observable permutation can be expressed as a combination of Hamiltonian permutations
            bool permutation_found = false;
            for (const auto& nullvector : nullspace) {
                if (nullvector.back() == 1) {  // Last element corresponds to observable permutation
                    permutation_found = true;
                    break;
                }
            }
            
            if (!permutation_found) {
                throw std::runtime_error("Error! The input observable cannot be written in terms of the permutation operators of the Hamiltonian");
            }
        }
        
        // Process observable groups
        std::vector<std::vector<bool>> obs_permutations;
        std::vector<std::pair<std::vector<bool>, std::vector<std::pair<std::complex<double>, std::vector<bool>>>>> sorted_obs_groups;
        
        for (auto& [perm, group] : obs_groups) {
            sorted_obs_groups.emplace_back(perm, std::move(group));
        }
        
        std::sort(sorted_obs_groups.begin(), sorted_obs_groups.end(),
                  [](const auto& a, const auto& b) {
                      return a.first < b.first;
                  });
        
        // Separate diagonal and off-diagonal terms for observable
        std::vector<std::complex<double>> obs_D0_coeff;
        std::vector<std::vector<bool>> obs_D0_product;
        std::vector<std::vector<bool>> obs_nontrivial_perms;
        std::vector<std::pair<std::vector<bool>, std::vector<std::pair<std::complex<double>, std::vector<bool>>>>> obs_offdiag_groups;
        
        for (const auto& [perm, group] : sorted_obs_groups) {
            bool is_identity = std::all_of(perm.begin(), perm.end(), [](bool x) { return !x; });
            
            if (is_identity) {
                // Handle diagonal terms
                for (const auto& [coeff, z_prod] : group) {
                    obs_D0_coeff.push_back(coeff);
                    obs_D0_product.push_back(z_prod);
                }
            } else {
                obs_nontrivial_perms.push_back(perm);
                obs_offdiag_groups.emplace_back(perm, group);
            }
        }
        
        // Store observable diagonal terms
        obs_data.MD0_coeff.push_back(obs_D0_coeff);
        obs_data.MD0_product.push_back(obs_D0_product);
        obs_data.MD0_size.back() = obs_D0_coeff.size();
        
        // Store observable off-diagonal terms
        std::vector<std::vector<std::complex<double>>> obs_D_coeff;
        std::vector<std::vector<std::vector<bool>>> obs_D_product;
        std::vector<int> obs_D_size;
        int obs_D_maxsize = 0;
        
        for (const auto& [perm, group] : obs_offdiag_groups) {
            std::vector<std::complex<double>> coeffs;
            std::vector<std::vector<bool>> products;
            
            for (const auto& [coeff, z_prod] : group) {
                coeffs.push_back(coeff);
                products.push_back(z_prod);
            }
            
            obs_D_coeff.push_back(coeffs);
            obs_D_product.push_back(products);
            obs_D_size.push_back(coeffs.size());
            obs_D_maxsize = std::max(obs_D_maxsize, (int)coeffs.size());
        }
        
        obs_data.MD_coeff.push_back(obs_D_coeff);
        obs_data.MD_product.push_back(obs_D_product);
        obs_data.MD_size.push_back(obs_D_size);
        obs_data.MD_maxsize.back() = obs_D_maxsize;
        
        // Create mapping matrix (MP) using nullspace method as in legacy prepare.cpp
        std::vector<std::vector<bool>> mapping_matrix;
        for (const auto& obs_perm : obs_nontrivial_perms) {
            // Convert Hamiltonian permutations to integer matrix
            std::vector<std::vector<int>> ham_perm_int;
            for (int i = 0; i < result.Nop; ++i) {
                std::vector<int> perm_row;
                for (bool bit : result.P_matrix[i]) {
                    perm_row.push_back(bit ? 1 : 0);
                }
                ham_perm_int.push_back(perm_row);
            }
            
            // Convert observable permutation to integer vector
            std::vector<int> obs_perm_int;
            for (bool bit : obs_perm) {
                obs_perm_int.push_back(bit ? 1 : 0);
            }
            
            // Create combined matrix for nullspace computation
            std::vector<std::vector<int>> combined_matrix = ham_perm_int;
            combined_matrix.push_back(obs_perm_int);
            
            // Compute nullspace
            auto nullspace = compute_nullspace(combined_matrix);
            
            // Find the nullspace vector with last element = 1 and minimum number of 1s
            std::vector<int> best_mapping(result.Nop, 0);
            bool found_mapping = false;
            int min_ones = result.Nop + 1;
            
            for (const auto& nullvector : nullspace) {
                if (nullvector.back() == 1) {
                    int num_ones = 0;
                    for (int i = 0; i < result.Nop; ++i) {
                        if (nullvector[i] == 1) num_ones++;
                    }
                    if (num_ones < min_ones) {
                        min_ones = num_ones;
                        best_mapping = nullvector;
                        found_mapping = true;
                    }
                }
            }
            
            // Convert to boolean mapping row (excluding the last element which corresponds to observable)
            std::vector<bool> mapping_row(result.Nop, false);
            if (found_mapping) {
                for (int i = 0; i < result.Nop; ++i) {
                    mapping_row[i] = (best_mapping[i] == 1);
                }
            }
            
            mapping_matrix.push_back(mapping_row);
        }
        
        obs_data.MP.push_back(mapping_matrix);
        obs_data.MNop.back() = obs_nontrivial_perms.size();
    }
    
    // Store observable data in result
    result.observable_data = std::move(obs_data);
    
    return result;
}

} // namespace pmrqmc