#include <OpSum.h>
#include <PMR.h>
#include <iostream>
#include <iomanip>
#include <cmath>

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
        for (int j = result.Nop - 1; j >= 0; --j) {
            std::cout << result.cycles[i][j];
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
    // Heisenberg model on 6 qubits
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

    // Observable 1: Simple Observable - Total Magnetization
    std::cout << "\n\n=== Observable 1: Simple Observable - Total Magnetization ===" << std::endl;
    
    OpSum mag("magnetization");
    for (int i = 0; i < N_QUBITS; ++i) {
        mag.add(1.0, "Z", i);
    }
    
    std::cout << "Magnetization observable:" << std::endl;
    mag.print();
    std::cout << std::endl;

    // Generate PMR for observable
    auto mag_pmr = pmr_obs(mag, ham_pmr);
    
    if (mag_pmr.has_observable_data()) {
        print_observable_data(mag_pmr.get_observable_data());
    }

    // Observable 2: Momentum Observable - Magnetic Structure Factor at q=π
    std::cout << "\n\n=== Observable 2: Momentum Observable - <|S^z(π)|²> ===" << std::endl;

    OpSum mom_obs("momentum_q_pi"); // S_z(π) S_z(-π), but since q=π ≡ -π mod 2π, it's |S^z(π)|²

    // S^z(q) = (1/√N) ∑_j e^{-i q r_j} S^z_j
    // For q = π, r_j = j (1D chain, unit spacing), e^{-i π j} = (-1)^j (real)
    // Thus <|S^z(π)|²> = <S^z(π) S^z(-π)> = (1/N) ∑_{j,k} e^{-i π (j - k)} <S^z_j S^z_k> = (1/N) ∑_{j,k} (-1)^{j-k} <S^z_j S^z_k>
    // Note: (-1)^{j-k} = (-1)^{j+k} (since (-1)^{-k} = (-1)^k), so phase based on j+k parity is equivalent
    const double norm_factor = 1.0 / static_cast<double>(N_QUBITS);

    // Naive O(N²) construction: Fine for small N_QUBITS; for large N, consider measuring single-body Fourier first, then squaring
    for (int j = 0; j < N_QUBITS; ++j) {
        for (int k = 0; k < N_QUBITS; ++k) {
            // Phase from (-1)^{j+k} ≡ (-1)^{j-k}
            int phase_parity = (j + k) % 2;  // 0: even (sign +1), 1: odd (sign -1)
            double sign = (phase_parity == 0) ? 1.0 : -1.0;
            double coeff = norm_factor * sign;
            // If S^z = Z/2, multiply coeff by 1/4 here; otherwise, assuming Z represents S^z directly
            mom_obs.add(coeff, "Z", j, "Z", k);
        }
    }

    std::cout << "Momentum observable S_z(π)S_z:" << std::endl;
    mom_obs.print();
    std::cout << std::endl;

    // Generate PMR for momentum observable
    auto mom_pmr = pmr_obs(mom_obs, ham_pmr);

    if (mom_pmr.has_observable_data()) {
        print_observable_data(mom_pmr.get_observable_data());
    }
    
    return 0;
}