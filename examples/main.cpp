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

int main() {
    // Example 1: Heisenberg model on 4 qubits
    std::cout << "=== Example 1: Heisenberg Model ===" << std::endl;
    
    OpSum hamiltonian("heisenberg");
    const int N_QUBITS = 4;
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
    
    // Example 2: Simple 2-qubit Hamiltonian
    std::cout << "\n\n=== Example 2: Simple 6c-Qubit Hamiltonian ===" << std::endl;
    
    OpSum ham2;
    ham2.add(1.000000, "Z", 0, "Z", 1);
    ham2.add(1.000000, "Z", 0, "Z", 2);
    ham2.add(1.000000, "Z", 1, "Z", 2);
    ham2.add(1.000000, "Z", 1, "Z", 3);
    ham2.add(1.000000, "Z", 2, "Z", 3);
    ham2.add(0.500000, "Z", 2, "Z", 4);
    ham2.add(0.500000, "Z", 1, "Z", 4);
    ham2.add(1.000000, "Z", 3, "Z", 4);
    ham2.add(1.000000, "Z", 3, "Z", 5);
    ham2.add(1.000000, "Z", 4, "Z", 5);
    ham2.add(1.000000, "Z", 2);
    ham2.add(-1.000000, "Z", 1);
    ham2.add(0.707107, "X", 0, "X", 1);
    ham2.add(0.707107, "Y", 0, "Y", 1);
    ham2.add(0.707107, "X", 0, "X", 1, "Z", 2);
    ham2.add(0.707107, "Y", 0, "Y", 1, "Z", 2);
    ham2.add(-0.707107, "X", 0, "X", 2, "Z", 1);
    ham2.add(-0.707107, "Y", 0, "Y", 2, "Z", 1);
    ham2.add(0.707107, "X", 0, "X", 2);
    ham2.add(0.707107, "Y", 0, "Y", 2);
    ham2.add(-0.500000, "X", 1, "X", 2, "Z", 4);
    ham2.add(-0.500000, "Y", 1, "Y", 2, "Z", 4);
    ham2.add(0.707107, "X", 1, "X", 3);
    ham2.add(0.707107, "Y", 1, "Y", 3);
    ham2.add(0.707107, "X", 1, "X", 3, "Z", 2);
    ham2.add(1.414214, "Y", 1, "Y", 3, "Z", 2);
    ham2.add(-0.707107, "X", 2, "X", 3, "Z", 1);
    ham2.add(-0.707107, "Y", 2, "Y", 3, "Z", 1);
    ham2.add(0.707107, "X", 2, "X", 3);
    ham2.add(1.414214, "Y", 2, "Y", 3);
    ham2.add(0.707107, "X", 1, "X", 4, "Z", 2);
    ham2.add(0.707107, "X", 2, "X", 4);
    ham2.add(1.000000, "X", 3, "X", 4);
    ham2.add(1.000000, "Y", 3, "Y", 4);
    ham2.add(1.000000, "X", 3, "X", 5);
    ham2.add(1.000000, "Y", 3, "Y", 5);
    ham2.add(1.000000, "X", 4, "X", 5);
    ham2.add(1.000000, "Y", 4, "Y", 5);

    ham2.print();
    std::cout << std::endl;
    
    auto ham2_pmr = pmr(ham2);
    print_pmr_result(ham2_pmr);

    // Example 3: Spin operators demonstration
    std::cout << "\n\n=== Example 3: Spin Operators ===" << std::endl;
    
    OpSum spin_ham("spin_demo");
    
    // Heisenberg model using spin operators (automatically handles 1/2 factor)
    const double J = 1.0;
    spin_ham.add(J, "Sx", 0, "Sx", 1);  // J * Sx_0 * Sx_1 = J * (X_0/2) * (X_1/2) = J/4 * X_0 X_1
    spin_ham.add(J, "Sy", 0, "Sy", 1);  // J * Sy_0 * Sy_1 = J * (Y_0/2) * (Y_1/2) = J/4 * Y_0 Y_1
    spin_ham.add(J, "Sz", 0, "Sz", 1);  // J * Sz_0 * Sz_1 = J * (Z_0/2) * (Z_1/2) = J/4 * Z_0 Z_1
    
    // Single spin operators
    spin_ham.add(2.0, "Sx", 2);        // 2.0 * Sx_2 = 2.0 * (X_2/2) = 1.0 * X_2
    spin_ham.add(1.0, "Sy", 2);        // 1.0 * Sy_2 = 1.0 * (Y_2/2) = 0.5 * Y_2
    spin_ham.add(3.0, "Sz", 2);        // 3.0 * Sz_2 = 3.0 * (Z_2/2) = 1.5 * Z_2
    
    // Only Sx, Sy, Sz operators supported
    // spin_ham.add(1.0, "S+", 0, "S-", 1); // S+ and S- removed for simplicity
    
    spin_ham.print();
    std::cout << std::endl;
    
    auto spin_pmr = pmr(spin_ham);
    print_pmr_result(spin_pmr);
    
    return 0;
}