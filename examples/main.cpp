#include <OpSum.h>
#include <PMR.h>
#include <iostream>
#include <iomanip>

using namespace pmrqmc;

void print_pmr_result(const PMR& result) {
    std::cout << "=== PMR Result ===" << std::endl;
    std::cout << "N (qubits): " << result.N << std::endl;
    std::cout << "Nop (permutation operators): " << result.Nop << std::endl;
    std::cout << "Ncycles (fundamental cycles): " << result.Ncycles << std::endl;
    std::cout << "D0_size (diagonal terms): " << result.D0_size << std::endl;
    std::cout << "D_maxsize (max off-diagonal terms): " << result.D_maxsize << std::endl;
    
    std::cout << "\n=== Permutation Matrices ===" << std::endl;
    for (int i = 0; i < result.Nop; ++i) {
        std::cout << "P[" << i << "]: ";
        for (bool bit : result.P_matrix[i]) {
            std::cout << bit;
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== Cycles ===" << std::endl;
    for (int i = 0; i < result.Ncycles; ++i) {
        std::cout << "Cycle[" << i << "]: ";
        for (bool bit : result.cycles[i]) {
            std::cout << bit;
        }
        std::cout << std::endl;
    }
    
    if (result.has_diagonal_terms()) {
        std::cout << "\n=== Diagonal Terms (D0) ===" << std::endl;
        for (int i = 0; i < result.D0_size; ++i) {
            std::cout << "D0_coeff[" << i << "]: " << result.D0_coeff[i] << std::endl;
            std::cout << "D0_product[" << i << "]: ";
            for (bool bit : result.D0_product[i]) {
                std::cout << bit;
            }
            std::cout << std::endl;
        }
    }
    
    if (result.has_offdiagonal_terms()) {
        std::cout << "\n=== Off-Diagonal Terms ===" << std::endl;
        for (int i = 0; i < result.Nop; ++i) {
            std::cout << "Operator " << i << " (" << result.D_size[i] << " terms):" << std::endl;
            for (int j = 0; j < result.D_size[i]; ++j) {
                std::cout << "  coeff[" << j << "]: " << result.D_coeff[i][j] << std::endl;
                std::cout << "  product[" << j << "]: ";
                for (bool bit : result.D_product[i][j]) {
                    std::cout << bit;
                }
                std::cout << std::endl;
            }
        }
    }
}

void print_observable_data(const ObservablePMR& obs_data) {
    std::cout << "=== Observable PMR Data ===" << std::endl;
    std::cout << "Number of observables: " << obs_data.get_num_observables() << std::endl;
    
    for (int obs_idx = 0; obs_idx < obs_data.get_num_observables(); ++obs_idx) {
        std::cout << "\n--- Observable " << obs_idx << ": " << obs_data.Mnames[obs_idx] << " ---" << std::endl;
        std::cout << "MNop[" << obs_idx << "]: " << obs_data.MNop[obs_idx] << std::endl;
        std::cout << "MD0_size[" << obs_idx << "]: " << obs_data.MD0_size[obs_idx] << std::endl;
        std::cout << "MD_maxsize[" << obs_idx << "]: " << obs_data.MD_maxsize[obs_idx] << std::endl;
        
        if (obs_data.has_diagonal_terms(obs_idx)) {
            std::cout << "\n  Diagonal Terms (MD0):" << std::endl;
            const auto& d0_coeffs = obs_data.get_diagonal_coeffs(obs_idx);
            for (size_t i = 0; i < d0_coeffs.size(); ++i) {
                std::cout << "    MD0_coeff[" << i << "]: " << d0_coeffs[i] << std::endl;
                std::cout << "    MD0_product[" << i << "]: ";
                for (bool bit : obs_data.MD0_product[obs_idx][i]) {
                    std::cout << bit;
                }
                std::cout << std::endl;
            }
        }
        
        if (obs_data.has_offdiagonal_terms(obs_idx)) {
            std::cout << "\n  Off-Diagonal Terms (MD):" << std::endl;
            const auto& md_coeffs = obs_data.get_offdiagonal_coeffs(obs_idx);
            for (size_t i = 0; i < md_coeffs.size(); ++i) {
                std::cout << "    MD_coeff[" << i << "] (" << obs_data.MD_size[obs_idx][i] << " terms):" << std::endl;
                for (size_t j = 0; j < obs_data.MD_size[obs_idx][i]; ++j) {
                    std::cout << "      coeff[" << j << "]: " << md_coeffs[i][j] << std::endl;
                    std::cout << "      product[" << j << "]: ";
                    for (bool bit : obs_data.MD_product[obs_idx][i][j]) {
                        std::cout << bit;
                    }
                    std::cout << std::endl;
                }
            }
        }
        
        // Show mapping matrix (MP)
        if (obs_data.MNop[obs_idx] > 0) {
            std::cout << "\n  Mapping Matrix (MP):" << std::endl;
            for (size_t i = 0; i < obs_data.MP[obs_idx].size(); ++i) {
                std::cout << "    MP[" << i << "]: ";
                for (bool bit : obs_data.MP[obs_idx][i]) {
                    std::cout << bit;
                }
                std::cout << std::endl;
            }
        }
    }
}


int main() {
    // Example 1: Heisenberg model on 4 qubits
    std::cout << "=== Example : Heisenberg Model ===" << std::endl;
    
    OpSum hamiltonian("heisenberg");
    const int N_QUBITS = 6;
    const double J_coupling = 1.0;
    const double h_field = 0.5;
    
    // Add Heisenberg interaction terms: -J * S_i . S_{i+1}
    for (int i = 0; i < N_QUBITS; ++i) {
        int next_i = (i + 1) % N_QUBITS;
        hamiltonian.add(-J_coupling, "X", i, "X", next_i);
        hamiltonian.add(-J_coupling, "Y", i, "Y", next_i);
        hamiltonian.add(-J_coupling, "Z", i, "Z", next_i);
    }
    
    // Add transverse field: -h * Σ X_i
    for (int i = 0; i < N_QUBITS; ++i) {
        hamiltonian.add(-h_field, "X", i);
    }
    
    hamiltonian.print();
    std::cout << std::endl;
    
    // Generate PMR representation
    auto ham_pmr = pmr(hamiltonian);
    print_pmr_result(ham_pmr);

    // Example 1: Single Observable - Magnetization
    std::cout << "\n\n=== Observable 1: Single Observable - Magnetization ===" << std::endl;
    
    OpSum mag("obs");
    for (int i = 0; i < N_QUBITS; ++i) {
        mag.add(1.0, "Z", i);
    }
    
    std::cout << "Magnetization observable:" << std::endl;
    mag.print();
    std::cout << std::endl;

    // Generate PMR for observable
    auto mag_pmr = ObservableBuilder::create_single(mag, ham_pmr);
    
    if (mag_pmr.has_observable_data()) {
        print_observable_data(mag_pmr.get_observable_data());
    }

    // Example 2: Spin-Spin Correlation Functions using OpSumBulk
    std::cout << "\n\n=== Example 2: Spin-Spin Correlation Functions ===" << std::endl;
    
    // Create bulk observables for correlation functions O_r = 1/N \sum_{r = r_i-r_j} Z_i Z_j
    OpSumBulk correlation_observables("spin_correlations");
    
    // Generate correlation observables for different distances r
    const int max_distance = N_QUBITS / 2;  // Maximum distance to consider
    
    for (int r = 1; r <= max_distance; ++r) {
        OpSum correlation_r("corr_r_" + std::to_string(r));
        
        // Sum over all pairs of sites with distance r (considering periodic boundary)
        for (int i = 0; i < N_QUBITS; ++i) {
            int j = (i + r) % N_QUBITS;  // Periodic boundary conditions
            correlation_r.add(1.0 / N_QUBITS, "Z", i, "Z", j);
        }// For now, we'll add normalization in the PMR processing
        
        correlation_observables.add(correlation_r);
    }
    
    // Generate PMR for all correlation observables
    auto correlations_pmr = ObservableBuilder::create_bulk(correlation_observables, ham_pmr);
    
    if (correlations_pmr.has_observable_data()) {
        print_observable_data(correlations_pmr.get_observable_data());
    }

    return 0;
}