#include <OpSum.h>
#include <iostream>
#include <iomanip>

using namespace pmrqmc;

void print_pmr_result(const PMRResult& result) {
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

void print_observable_data(const ObservablePMRResult& obs_data) {
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
    // Create the same Hamiltonian as Example 2 in main.cpp
    std::cout << "=== Creating Hamiltonian (from main.cpp Example 2) ===" << std::endl;
    
    OpSum hamiltonian("legacy_hamiltonian");
    
    // Add all the terms from the legacy Hamiltonian
    hamiltonian.add(1.000000, "Z", 0, "Z", 1);
    hamiltonian.add(1.000000, "Z", 0, "Z", 2);
    hamiltonian.add(1.000000, "Z", 1, "Z", 2);
    hamiltonian.add(1.000000, "Z", 1, "Z", 3);
    hamiltonian.add(1.000000, "Z", 2, "Z", 3);
    hamiltonian.add(0.500000, "Z", 2, "Z", 4);
    hamiltonian.add(0.500000, "Z", 1, "Z", 4);
    hamiltonian.add(1.000000, "Z", 3, "Z", 4);
    hamiltonian.add(1.000000, "Z", 3, "Z", 5);
    hamiltonian.add(1.000000, "Z", 4, "Z", 5);
    hamiltonian.add(1.000000, "Z", 2);
    hamiltonian.add(-1.000000, "Z", 1);
    hamiltonian.add(0.707107, "X", 0, "X", 1);
    hamiltonian.add(0.707107, "Y", 0, "Y", 1);
    hamiltonian.add(0.707107, "X", 0, "X", 1, "Z", 2);
    hamiltonian.add(0.707107, "Y", 0, "Y", 1, "Z", 2);
    hamiltonian.add(-0.707107, "X", 0, "X", 2, "Z", 1);
    hamiltonian.add(-0.707107, "Y", 0, "Y", 2, "Z", 1);
    hamiltonian.add(0.707107, "X", 0, "X", 2);
    hamiltonian.add(0.707107, "Y", 0, "Y", 2);
    hamiltonian.add(-0.500000, "X", 1, "X", 2, "Z", 4);
    hamiltonian.add(-0.500000, "Y", 1, "Y", 2, "Z", 4);
    hamiltonian.add(0.707107, "X", 1, "X", 3);
    hamiltonian.add(0.707107, "Y", 1, "Y", 3);
    hamiltonian.add(0.707107, "X", 1, "X", 3, "Z", 2);
    hamiltonian.add(1.414214, "Y", 1, "Y", 3, "Z", 2);
    hamiltonian.add(-0.707107, "X", 2, "X", 3, "Z", 1);
    hamiltonian.add(-0.707107, "Y", 2, "Y", 3, "Z", 1);
    hamiltonian.add(0.707107, "X", 2, "X", 3);
    hamiltonian.add(1.414214, "Y", 2, "Y", 3);
    hamiltonian.add(0.707107, "X", 1, "X", 4, "Z", 2);
    hamiltonian.add(0.707107, "X", 2, "X", 4);
    hamiltonian.add(1.000000, "X", 3, "X", 4);
    hamiltonian.add(1.000000, "Y", 3, "Y", 4);
    hamiltonian.add(1.000000, "X", 3, "X", 5);
    hamiltonian.add(1.000000, "Y", 3, "Y", 5);
    hamiltonian.add(1.000000, "X", 4, "X", 5);
    hamiltonian.add(1.000000, "Y", 4, "Y", 5);

    std::cout << "Hamiltonian created with " << hamiltonian.get_max_qubit_index() << " qubits" << std::endl;
    std::cout << std::endl;
    
    // Generate PMR for Hamiltonian
    auto ham_pmr = pmr(hamiltonian);
    print_pmr_result(ham_pmr);
    
    // Example 1: Single Observable - Magnetization
    std::cout << "\n\n=== Example 1: Single Observable - Magnetization ===" << std::endl;
    
    OpSum obstest("obs");
    obstest.add(1.000000, "Z", 0, "Z", 1);
    obstest.add(1.000000, "Z", 0, "Z", 2);
    obstest.add(1.000000, "Z", 1, "Z", 2);
    obstest.add(1.000000, "Z", 1, "Z", 3);
    obstest.add(1.000000, "Z", 2, "Z", 3);
    obstest.add(0.500000, "Z", 2, "Z", 4);
    obstest.add(0.500000, "Z", 1, "Z", 4);
    obstest.add(1.000000, "Z", 3, "Z", 4);
    obstest.add(1.000000, "Z", 3, "Z", 5);
    obstest.add(1.000000, "Z", 4, "Z", 5);
    obstest.add(1.000000, "Z", 2);
    obstest.add(-1.000000, "Z", 1);
    
    std::cout << "Magnetization observable:" << std::endl;
    obstest.print();
    std::cout << std::endl;
    
    // Generate PMR for observable
    auto mag_pmr = ObservableBuilder::create_single(obstest, ham_pmr);
    
    if (mag_pmr.has_observable_data()) {
        print_observable_data(mag_pmr.get_observable_data());
    }

    // Example 2: Single Observable - Magnetization
    std::cout << "\n\n=== Example 1: Single Observable - Magnetization ===" << std::endl;
    
    OpSum obstest2("obs");
    obstest2.add(0.707107, "X", 0, "X", 1);
    obstest2.add(0.707107, "Y", 0, "Y", 1);
    obstest2.add(0.707107, "X", 0, "X", 1, "Z", 2);
    obstest2.add(0.707107, "Y", 0, "Y", 1, "Z", 2);
    obstest2.add(-0.707107, "X", 0, "X", 2, "Z", 1);
    obstest2.add(-0.707107, "Y", 0, "Y", 2, "Z", 1);
    obstest2.add(0.707107, "X", 0, "X", 2);
    obstest2.add(0.707107, "Y", 0, "Y", 2);
    obstest2.add(-0.500000, "X", 1, "X", 2, "Z", 4);
    obstest2.add(-0.500000, "Y", 1, "Y", 2, "Z", 4);
    obstest2.add(0.707107, "X", 1, "X", 3);
    obstest2.add(0.707107, "Y", 1, "Y", 3);
    obstest2.add(0.707107, "X", 1, "X", 3, "Z", 2);
    obstest2.add(1.414214, "Y", 1, "Y", 3, "Z", 2);
    obstest2.add(-0.707107, "X", 2, "X", 3, "Z", 1);
    obstest2.add(-0.707107, "Y", 2, "Y", 3, "Z", 1);
    obstest2.add(0.707107, "X", 2, "X", 3);
    obstest2.add(1.414214, "Y", 2, "Y", 3);
    obstest2.add(0.707107, "X", 1, "X", 4, "Z", 2);
    obstest2.add(0.707107, "X", 2, "X", 4);
    obstest2.add(1.000000, "X", 3, "X", 4);
    obstest2.add(1.000000, "Y", 3, "Y", 4);
    obstest2.add(1.000000, "X", 3, "X", 5);
    obstest2.add(1.000000, "Y", 3, "Y", 5);
    obstest2.add(1.000000, "X", 4, "X", 5);
    obstest2.add(1.000000, "Y", 4, "Y", 5);
    
    std::cout << "Magnetization observable:" << std::endl;
    obstest2.print();
    std::cout << std::endl;
    
    // Generate PMR for observable
    auto mag_pmr2 = ObservableBuilder::create_single(obstest2, ham_pmr);
    
    if (mag_pmr2.has_observable_data()) {
        print_observable_data(mag_pmr2.get_observable_data());
    }
    
    return 0;
}