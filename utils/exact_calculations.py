"""
This code was written in support of the experiments carried out in:
* Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
* Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).

Description: Contains code to do exact numerical computations with which to compare QMC outputs.
"""
# %%
################################################
# Functions to compute exact values for PRL
# model: [10.1103/PhysRevLett.100.100501]
# This is Mathematica code converted
# to Python using chatGPT4-o.
################################################

import numpy as np
import scipy

def get_paulis():
    """
    Gets 1 qubit Pauli basis.
    """
    id2 = np.array([[1, 0], [0, 1]], dtype=complex)
    px = np.array([[0, 1], [1, 0]], dtype=complex)
    py = np.array([[0, -1j], [1j, 0]], dtype=complex)
    pz = np.array([[1, 0], [0, -1]], dtype=complex)

    return id2, px, py, pz

def prl_h(Bz):
    """
    Returns 2 qubit PRL Hamiltonian.
    """
    id2, px, _, pz = get_paulis()
    h0 = np.kron(pz, pz) + 0.1 * (np.kron(px, id2) + np.kron(id2, px))
    h1 = np.kron(pz, id2) + np.kron(id2, pz)
    h = h0 + Bz * h1

    return h0, h1, h

def prl_gs_gtau(Bz, tau):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = <psi(Bz)|H_1|psi(Bz)> for |psi(Bz)>
    the ground-state of PRL model.
    """
    # Compute the eigenvalues and eigenvectors of H(Bz)
    _, h1, h = prl_h(Bz)
    evals, evecs = np.linalg.eig(h)
    
    # Find the ground state energy and corresponding eigenvector
    e0 = np.min(evals)
    psi0 = evecs[:, np.argmin(evals)]
    
    # Compute the correlation term
    corr = np.sum([
        np.exp(tau * (e0 - evals[j])) * np.abs(np.dot(np.conj(psi0.T), np.dot(h1, evecs[:, j])))**2
        for j in range(len(evals))
    ])
    
    # Compute the average value of H1 in the ground state
    avgH1 = np.dot(np.conj(psi0.T), np.dot(h1, psi0))
    
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def prl_gs_chiE(Bz):
    """
    Integrates prl_gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: prl_gs_gtau(Bz, tau), 0, np.inf)
    return result, error

def prl_gs_chiF(Bz):
    """
    Integrates tau*prl_gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * prl_gs_gtau(Bz, tau), 0, np.inf)
    return result, error

def prl_beta_gtau(Bz, beta, tau):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = Tr[H_1 rho(Bz, beta)] for rho the
    thermal state PRL model at inverse temp beta.
    """
    # Compute the eigenvalues and eigenvectors of H(Bz)
    _, h1, h = prl_h(Bz)
    evals, evecs = np.linalg.eig(h)
    
    # Compute the partition function Z
    z = np.sum(np.exp(-beta * evals))
    
    # Compute the correlation term
    corr = np.sum([
        np.exp(- (beta - tau) * evals[i]) * np.exp(- tau * evals[j]) *
        np.abs(np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, j])))**2
        for i in range(len(evals))
        for j in range(len(evals))
    ]) / z
    
    # Compute the average value of H1
    avgH1 = np.sum([
        np.exp(-beta * evals[i]) * np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, i]))
        for i in range(len(evals))
    ]) / z
    
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def prl_beta_chiE(Bz, beta):
    """
    Integrates prl_beta_gtau over tau from [0, beta].
    """
    result, error = scipy.integrate.quad(lambda tau: prl_beta_gtau(Bz, beta, tau), 0, beta)
    return result, error

def prl_beta_chiX(Bz, beta):
    """
    Integrates tau*prl_beta_gtau over tau from [0, beta].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * prl_beta_gtau(Bz, beta, tau), 0, beta)
    return result, error

def prl_beta_chiF(Bz, beta):
    """
    Integrates tau*prl_beta_gtau over tau from [0, beta/2].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * prl_beta_gtau(Bz, beta, tau), 0, beta/2)
    return result, error

def prl_gs_fidsus(Bz, epsilon=0.0001):
    """
    Returns T = 0 fid sus for PRL model.
    """
    # Compute the ground state vector for H(Bz)
    _, _, h = prl_h(Bz)
    _, evecs = scipy.sparse.linalg.eigsh(h, k=1, which='SA')
    psi0 = evecs[:, 0]
    
    # Compute the ground state vector for H(Bz + epsilon)
    _, _, heps = prl_h(Bz + epsilon)
    _, evecs_eps = scipy.sparse.linalg.eigsh(heps, k=1, which='SA')
    psiPeps = evecs_eps[:, 0]
    
    # Compute dPsi
    dPsi = (psiPeps - psi0) / epsilon
    
    # Compute chi
    chi = (np.dot(np.conj(dPsi.T), dPsi) - 
           np.dot(np.conj(dPsi.T), psi0) * np.dot(np.conj(psi0.T), dPsi))
    if chi.imag < 1e-8:
        chi = chi.real

    return chi
# %%
# %%

# %%
########################################################
# Exact computation for H.txt: D_0 = H_1 case
# This code was partially generated with chatGPT4-0
# and verified carefully by hand and tests.
########################################################
import numpy as np
import scipy.sparse as sp


#--------------------------------------------
# Some I/O helpers
#--------------------------------------------
def sparse_pauli(pauli):
    """
    Return a sparse matrix representation of the given Pauli operator.
    
    Args:
        pauli (str): Pauli operator ('I', 'X', 'Y', 'Z').
    
    Returns:
        scipy.sparse.csr_matrix: Sparse matrix representation of the Pauli operator.
    """
    if pauli == 'I':
        return sp.csr_matrix(np.eye(2))
    elif pauli == 'X':
        return sp.csr_matrix([[0, 1], [1, 0]])
    elif pauli == 'Y':
        return sp.csr_matrix([[0, -1j], [1j, 0]])
    elif pauli == 'Z':
        return sp.csr_matrix([[1, 0], [0, -1]])
    else:
        raise ValueError(f"Unknown Pauli operator: {pauli}")

#def load_hamiltonian(num_qubits, filename="/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/tfim_experiments/H.txt"):
#def load_hamiltonian(num_qubits, filename="../H.txt"):
def load_hamiltonian(num_qubits, filename="/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/H.txt"):
    """
    Load a Hamiltonian from a file and return it as a sparse matrix.
    
    Args:
        filename (str): Path to the file describing the Hamiltonian.
        num_qubits (int): Number of qubits in the system.
        
    Returns:
        scipy.sparse.csr_matrix: The Hamiltonian as a sparse matrix.
    """
    # Initialize Hamiltonian as a sparse zero matrix
    hamiltonian = sp.csr_matrix((2**num_qubits, 2**num_qubits))

    with open(filename, 'r') as file:
        for line in file:
            # Parse the line
            parts = line.strip().split()
            coefficient = float(parts[0])  # Extract the coefficient
            operators = parts[1:]          # Extract the operators
            
            # Start with the identity matrix for the full system
            term = sp.eye(1, format="csr")
            
            for i in range(num_qubits):
                # Check if a specific qubit has a Pauli operator; default to identity
                op = 'I'
                for j in range(0, len(operators), 2):
                    if int(operators[j]) - 1 == i:  # Adjust qubit index from 1-based to 0-based
                        op = operators[j+1]
                        break

                # Create the operator for the current qubit
                single_qubit_op = sparse_pauli(op)
                
                # Expand to the full system
                term = sp.kron(term, single_qubit_op, format="csr")

            # Add the term to the Hamiltonian
            hamiltonian += coefficient * term

    return hamiltonian

#--------------------------------------------
# Exact groundstate calculations
#--------------------------------------------
def htxt_diagonal_exact_gs_gtau(num_qubits, tau):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = <psi(Bz)|H_1|psi(Bz)> for |psi(Bz)>
    the ground-state of H.txt Hamiltonian.
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    evals, evecs = scipy.linalg.eigh(h.toarray())

    # Find the ground state energy and corresponding eigenvector
    e0 = np.min(evals)
    psi0 = evecs[:, np.argmin(evals)]

    # Compute the average value of H1 in the ground state
    avgH1 = psi0.conj().transpose() @ h1 @ psi0
    
    # Compute the correlation term
    corr = 0.0
    for j in range(len(evals)):
        coeff = np.exp(tau * (e0 - evals[j]))
        h1_vec = h1 @ evecs[:, j]
        exp_val = np.abs(psi0.conj().transpose() @ h1_vec)
        corr += coeff * exp_val
    
    gTau = corr - (avgH1)**2
    # Truncate small imaginary part
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def htxt_diagonal_exact_gs_chiE(n):
    """
    Integrates exact, diagonal gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: htxt_diagonal_exact_gs_gtau(n, tau), 0, np.inf)
    return result, error

def htxt_diagonal_exact_gs_chiF(n):
    """
    Integrates exact, diagonal tau*gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * htxt_diagonal_exact_gs_gtau(n, tau), 0, np.inf)
    return result, error

#--------------------------------------------
# Exact thermal state calculations
#--------------------------------------------
def htxt_diagonal_exact_beta_h1(n, beta):
    """
    Returns <H_1> for <H_1> = Tr[H_1 rho(Bz, beta)] 
    and rho the thermal state of Htxt at inverse
    temperature \beta.
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(n)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr').toarray()
    evals, evecs = scipy.linalg.eigh(h.toarray())
    
    # Compute the partition function Z
    z = np.sum(np.exp(-beta * evals))

    # Compute the average value of H1
    avgH1 = np.sum([
        np.exp(-beta * evals[i]) * np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, i]))
        for i in range(len(evals))
    ]) / z
    # make real by truncating small imaginary parts
    if avgH1.imag < 1e-8:
        avgH1 = avgH1.real

    return avgH1

def htxt_diagonal_exact_beta_corr(n, beta, tau):
    """
    Returns <H_1(tau)H_1> for <H_1> = Tr[H_1 rho(Bz, beta)] 
    and rho the thermal state of Htxt at inverse
    temperature \beta.
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(n)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr').toarray()
    evals, evecs = scipy.linalg.eigh(h.toarray())
    
    # Compute the partition function Z
    z = np.sum(np.exp(-beta * evals))

    # Compute the correlation term
    corr = np.sum([
        np.exp(- (beta - tau) * evals[i]) * np.exp(- tau * evals[j]) *
        np.abs(np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, j])))**2
        for i in range(len(evals))
        for j in range(len(evals))
    ]) / z
    # make real by truncating small imaginary parts
    if corr.imag < 1e-8:
        corr = corr.real

    return corr

def htxt_diagonal_exact_beta_gtau(n, beta, tau):
    """
    Returns <H_1(tau)H_1> for <H_1> = Tr[H_1 rho(Bz, beta)] 
    and rho the thermal state of Htxt at inverse
    temperature \beta.
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(n)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr').toarray()
    evals, evecs = scipy.linalg.eigh(h.toarray())
    
    # Compute the partition function Z
    z = np.sum(np.exp(-beta * evals))

    # Compute the average value of H1
    avgH1 = np.sum([
        np.exp(-beta * evals[i]) * np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, i]))
        for i in range(len(evals))
    ]) / z

    # Compute the correlation term
    corr = np.sum([
        np.exp(- (beta - tau) * evals[i]) * np.exp(- tau * evals[j]) *
        np.abs(np.dot(np.conj(evecs[:, i].T), np.dot(h1, evecs[:, j])))**2
        for i in range(len(evals))
        for j in range(len(evals))
    ]) / z
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def htxt_diagonal_exact_beta_chiEint(n, beta):
    """
    Integrates htxt_diagonal_exact_beta_gtau over tau from [0, beta].
    """
    result, error = scipy.integrate.quad(lambda tau: htxt_diagonal_exact_beta_corr(n, beta, tau), 0, beta)
    return result, error

def htxt_diagonal_exact_beta_chiFint(n, beta):
    """
    Integrates tau*htxt_diagonal_exact_beta_gtau over tau from [0, beta/2].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * htxt_diagonal_exact_beta_corr(n, beta, tau), 0, beta/2)
    return result, error

def htxt_diagonal_exact_fs_experiment_obs(n, beta):
    """
    Returns the three output observables in a FS experiment run
    for the current Hamiltonian in H.txt; namely,
    - <H1>
    - \int_0^\beta <H1(\tau)H1> \dtau (ES correlation integral)
    - \int_0^{\beta/2} \tau <H1(\tau)H1> \dtau (FS correlation integral)
    """
    avgH1 = htxt_diagonal_exact_beta_h1(n, beta)
    intES = htxt_diagonal_exact_beta_chiEint(n, beta)
    intFS = htxt_diagonal_exact_beta_chiFint(n, beta)

    return avgH1, intES, intFS

#--------------------------------------------
# Truncated groundstate calculations
#--------------------------------------------
def htxt_diagonal_approx_gs_gtau(num_qubits, tau, trunc=5):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = <psi(Bz)|H_1|psi(Bz)> for |psi(Bz)>
    the ground-state of H.txt Hamiltonian.
    """
    # get get first [trunc] eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    evals, evecs = scipy.sparse.linalg.eigsh(h, k=trunc, which='SA')

    # Find the ground state energy and corresponding eigenvector
    e0 = np.min(evals)
    psi0 = evecs[:, np.argmin(evals)]
    
    # Compute the correlation term
    corr = 0.0
    for j in range(len(evals)):
        coeff = np.exp(tau * (e0 - evals[j]))
        h1_vec = h1 @ evecs[:, j]
        exp_val = np.abs(psi0.conj().transpose() @ h1_vec)
        corr += coeff * exp_val

    # Compute the average value of H1 in the ground state
    avgH1 = psi0.conj().transpose() @ h1 @ psi0
    
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def htxt_diagonal_approx_gs_chiE(n, trunc=5):
    """
    Integrates approx, diagonal gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: htxt_diagonal_approx_gs_gtau(n, tau, trunc), 0, np.inf)
    return result, error

def htxt_diagonal_approx_gs_chiF(n, trunc=5):
    """
    Integrates approx, diagonal tau*gs_gtau over tau from [0, np.inf].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * htxt_diagonal_approx_gs_gtau(n, tau, trunc), 0, np.inf)
    return result, error

def htxt_diagonal_exact_fidsus(n, eps=0.00001):
    """
    Returns T = 0 fid sus for htxt model.
    """
    # Compute the ground state vector for H(lam)
    h = load_hamiltonian(n)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    _, evecs = scipy.sparse.linalg.eigsh(h, k=1, which='SA')
    psi0 = evecs[:, 0]
    
    # Compute the ground state vector for H(lam + epsilon)
    eps_h = h + eps * h1
    _, evecs_eps = scipy.sparse.linalg.eigsh(eps_h, k=1, which='SA')
    psiPeps = evecs_eps[:, 0]
    
    # Compute dPsi
    dPsi = (psiPeps - psi0) / eps
    
    # Compute chi
    dPsi_dag = dPsi.conj().transpose()
    psi0_dag = psi0.conj().transpose()
    chi = dPsi_dag @ dPsi - dPsi_dag @ psi0 * psi0_dag @ dPsi
    #chi = (np.dot(np.conj(dPsi.T), dPsi) - 
    #       np.dot(np.conj(dPsi.T), psi0) * np.dot(np.conj(psi0.T), dPsi))
    if chi.imag < 1e-8:
        chi = chi.real

    return chi

def htxt_diagonal_exact_esus(n, eps=0.00001):
    """
    Returns T = 0 energy sus for htxt model.
    """
    # Compute the ground state energy for H(lam)
    h = load_hamiltonian(n)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    evals, _ = scipy.sparse.linalg.eigsh(h, k=1, which='SA')
    e0 = evals[0]

    # Compute the ground state energy for H(lam + epsilon)
    forward_h = h + eps * h1
    forward_evals, _ = scipy.sparse.linalg.eigsh(forward_h, k=1, which='SA')
    forward_e0 = forward_evals[0]

    # Compute the ground state energy for H(lam - epsilon)
    backward_h = h - eps * h1
    backward_evals, _ = scipy.sparse.linalg.eigsh(backward_h, k=1, which='SA')
    backward_e0 = backward_evals[0]
    esus = -(forward_e0 - 2 * e0 + backward_e0) / (eps**2)

    return esus


#--------------------------------------------
# Truncated thermal state calculations
#--------------------------------------------
def htxt_diagonal_approx_beta_gtau(num_qubits, beta, tau, trunc=5):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = <psi(Bz)|H_1|psi(Bz)> for |psi(Bz)>
    the ground-state of H.txt Hamiltonian.
    """
    # get get first [trunc] eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    evals, evecs = scipy.sparse.linalg.eigsh(h, k=trunc, which='SA')

    # Compute the partition function Z
    z = np.sum(np.exp(-beta * evals))

    # Compute the average value of H1
    avgH1 = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        h1_ele = evecs[:, i].conj().transpose() @ h1 @ evecs[:, i]
        avgH1 += expi * h1_ele 
    avgH1 /= z

    # Compute the correlation term
    corr = 0.0
    for i in range(len(evals)):
        expi = np.exp(- (beta - tau) * evals[i])
        eveci_dag = evecs[:, i].conj().transpose()
        for j in range(len(evals)):
            expj = np.exp(- tau * evals[j])
            h1_ele = np.abs(eveci_dag @ h1 @ evecs[:, j])**2
            corr += expi * expj * h1_ele
    corr /= z
    
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def htxt_diagonal_approx_beta_chiE(n, beta, trunc=5):
    """
    Integrates htxt_diagonal_approx_beta_gtau over tau from [0, beta].
    """
    result, error = scipy.integrate.quad(lambda tau: htxt_diagonal_approx_beta_gtau(n, beta, tau, trunc), 0, beta)
    return result, error

def htxt_diagonal_approx_beta_chiF(n, beta, trunc=5):
    """
    Integrates tau*htxt_diagonal_exact_beta_gtau over tau from [0, beta/2].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * htxt_diagonal_approx_beta_gtau(n, beta, tau, trunc), 0, beta/2)
    return result, error

def htxt_diagonal_approx_proj_beta_gtau(num_qubits, beta, tau, proj='+', trunc=5):
    """
    Returns G(tau) = <H_1(tau)H_1> - <H_1>^2 for
    <H_1> = <psi(Bz)|H_1|psi(Bz)> for |psi(Bz)>
    the ground-state of H.txt Hamiltonian.
    """
    # get get first [trunc] eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr')
    evals, evecs = scipy.sparse.linalg.eigsh(h, k=trunc, which='SA')

    # define projector into plus or minus space
    z = sparse_pauli('Z')
    tensor_z = sp.eye(1, format="csr")
    for j in range(num_qubits):
        tensor_z = sp.kron(tensor_z, z, format="csr")
    id = sp.eye(2**num_qubits, format="csr")
    if proj == '+':
        proj = (id + tensor_z) / 2
    else:
        proj = (id - tensor_z) / 2

    # Compute the partition function Z
    z = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        proj_val = evecs[:, i].conj().transpose() @ proj @ evecs[:, i]
        z += expi * proj_val
    z = z.real

    # Compute the average value of H1
    avgH1 = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        h1_ele = evecs[:, i].conj().transpose() @ h1 @ proj @ evecs[:, i]
        avgH1 += expi * h1_ele 
    avgH1 /= z

    # Compute the correlation term
    corr = 0.0
    for i in range(len(evals)):
        expi = np.exp(- (beta - tau) * evals[i])
        eveci_dag = evecs[:, i].conj().transpose()
        for j in range(len(evals)):
            expj = np.exp(- tau * evals[j])
            h1_ele = np.abs(eveci_dag @ proj @ h1 @ evecs[:, j])**2
            corr += expi * expj * h1_ele
    corr /= z
    
    # Compute gTau
    gTau = corr - (avgH1)**2
    if gTau.imag < 1e-8:
        gTau = gTau.real

    return gTau

def htxt_diagonal_approx_proj_beta_chiE(n, beta, proj="+", trunc=5):
    """
    Integrates htxt_diagonal_approx_beta_gtau over tau from [0, beta].
    """
    result, error = scipy.integrate.quad(lambda tau: htxt_diagonal_approx_proj_beta_gtau(n, beta, tau, proj, trunc), 0, beta)
    return result, error

def htxt_diagonal_approx_proj_beta_chiF(n, beta, proj="+", trunc=5):
    """
    Integrates tau*htxt_diagonal_exact_beta_gtau over tau from [0, beta/2].
    """
    result, error = scipy.integrate.quad(lambda tau: tau * htxt_diagonal_approx_proj_beta_gtau(n, beta, tau, proj, trunc), 0, beta/2)
    return result, error
# %%
# %%
def htxt_lanczos(num_qubits, trunc):
    # get get first [trunc] eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    evals, evecs = scipy.sparse.linalg.eigsh(h, k=trunc, which='SA')

    return evals, evecs

# %%

# %%
#--------------------------------------------
# Exact, projected thermal state calculations
#--------------------------------------------
def htxt_diagonal_exact_proj_fs_experiment_obs(num_qubits, beta, proj="+"):
    """
    Computes fidelity susceptibility experiment quantities
    * <H_1>
    * ES integral, \int_0^\beta <H_1(t)H_1> dt
    * FS integral, \int_0^{\beta/2} t <H_1(t)H_1> dt
    restricted to positive or negative parity subspace. 
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = scipy.sparse.diags(h.diagonal(), format='csr').toarray()
    evals, evecs = scipy.linalg.eigh(h.toarray())
    
    # define projector into plus or minus space
    z = sparse_pauli('Z')
    tensor_z = sp.eye(1, format="csr")
    for j in range(num_qubits):
        tensor_z = sp.kron(tensor_z, z, format="csr")
    id = sp.eye(2**num_qubits, format="csr")
    if proj == '+':
        proj = (id + tensor_z) / 2
    else:
        proj = (id - tensor_z) / 2

    # Compute the partition function Z (no projection)
    z = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        # note one can add up only Z in correct parity subspace OR renorm by <parity> at end
        #proj_val = evecs[:, i].conj().transpose() @ proj @ evecs[:, i]
        proj_val = evecs[:, i].conj().transpose() @ evecs[:, i]
        z += expi * proj_val
    z = z.real

    # compute average parity
    # note one can add up only Z in correct parity subspace OR renorm by <parity> at end
    parity = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        parity_ele = evecs[:, i].conj().transpose() @ proj @ evecs[:, i]
        parity += expi * parity_ele 
    parity /= z
    parity = parity.real

    # Compute the average value of H1
    avgH1 = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        h1_ele = evecs[:, i].conj().transpose() @ h1 @ proj @ evecs[:, i]
        avgH1 += expi * h1_ele 
    avgH1 /= z
    avgH1 = avgH1.real / parity

    def corr_func(tau):
        corr = 0.0
        for i in range(len(evals)):
            expi = np.exp(- (beta - tau) * evals[i])
            eveci_dag = evecs[:, i].conj().transpose()
            for j in range(len(evals)):
                expj = np.exp(- tau * evals[j])
                h1_ele = np.abs(eveci_dag @ proj @ h1 @ evecs[:, j])**2
                corr += expi * expj * h1_ele
        corr /= z
        return corr.real / parity
        return corr.real

    intES, errorES = scipy.integrate.quad(corr_func, 0, beta)
    intFS, errorFS = scipy.integrate.quad(lambda tau: tau * corr_func(tau), 0, beta/2)

    return avgH1, intES, intFS, parity

#--------------------------------------------
# Exact, projected thermal state calculations
#--------------------------------------------
def htxt_offdiagonal_exact_proj_fs_experiment_obs(num_qubits, beta, proj="+"):
    """
    Computes fidelity susceptibility experiment quantities
    * <H_1>
    * ES integral, \int_0^\beta <H_1(t)H_1> dt
    * FS integral, \int_0^{\beta/2} t <H_1(t)H_1> dt
    restricted to positive or negative parity subspace. 
    """
    # get eigenvalues and eigenvectors
    h = load_hamiltonian(num_qubits)
    h1 = h.toarray() - scipy.sparse.diags(h.diagonal(), format='csr').toarray()
    evals, evecs = scipy.linalg.eigh(h.toarray())
    
    # define projector into plus or minus space
    z = sparse_pauli('Z')
    tensor_z = sp.eye(1, format="csr")
    for j in range(num_qubits):
        tensor_z = sp.kron(tensor_z, z, format="csr")
    id = sp.eye(2**num_qubits, format="csr")
    if proj == '+':
        proj = (id + tensor_z) / 2
    else:
        proj = (id - tensor_z) / 2

    # Compute the partition function Z (no projection)
    z = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        # note one can add up only Z in correct parity subspace OR renorm by <parity> at end
        #proj_val = evecs[:, i].conj().transpose() @ proj @ evecs[:, i]
        proj_val = evecs[:, i].conj().transpose() @ evecs[:, i]
        z += expi * proj_val
    z = z.real

    # compute average parity
    # note one can add up only Z in correct parity subspace OR renorm by <parity> at end
    parity = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        parity_ele = evecs[:, i].conj().transpose() @ proj @ evecs[:, i]
        parity += expi * parity_ele 
    parity /= z
    parity = parity.real

    # Compute the average value of H1
    avgH1 = 0.0
    for i in range(len(evals)):
        expi = np.exp(-beta * evals[i])
        h1_ele = evecs[:, i].conj().transpose() @ h1 @ proj @ evecs[:, i]
        avgH1 += expi * h1_ele 
    avgH1 /= z
    avgH1 = avgH1.real / parity

    def corr_func(tau):
        corr = 0.0
        for i in range(len(evals)):
            expi = np.exp(- (beta - tau) * evals[i])
            eveci_dag = evecs[:, i].conj().transpose()
            for j in range(len(evals)):
                expj = np.exp(- tau * evals[j])
                h1_ele = np.abs(eveci_dag @ proj @ h1 @ evecs[:, j])**2
                corr += expi * expj * h1_ele
        corr /= z
        return corr.real / parity
        return corr.real

    intES, errorES = scipy.integrate.quad(corr_func, 0, beta)
    intFS, errorFS = scipy.integrate.quad(lambda tau: tau * corr_func(tau), 0, beta/2)

    return avgH1, intES, intFS, parity
# %%

# %%

def sparse_kron_n(mat_list):
    """Kronecker product for a list of matrices."""
    result = mat_list[0]
    for mat in mat_list[1:]:
        result = sp.kron(result, mat, format="csr")
    return result

def sparse_tilted_ising_field_ham(n, theta, B):
    """
    
    """
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I = np.eye(2, dtype=complex)

    # Initialize ZZ Hamiltonian terms as a sparse matrix
    Hzz = sp.csr_matrix((2**n, 2**n), dtype=complex)
    # Add ZZ terms: \sum_{i=1}^{n-1} Z_i Z_{i+1}
    for i in range(n - 1):
        operators = [I] * n
        operators[i] = Z
        operators[i + 1] = Z
        Hzz += sparse_kron_n(operators)
    
    # Initialize local Hamiltonian terms as a sparse matrix
    Hdriver = sp.csr_matrix((2**n, 2**n), dtype=complex)
    # Add B(X cos(theta) + Z sin(theta)) terms: \sum_{i=1}^n
    for i in range(n):
        # X_i cos(theta)
        operators = [I] * n
        operators[i] = X
        Hdriver += B * np.cos(theta) * sparse_kron_n(operators)
        
        # Z_i sin(theta)
        operators[i] = Z
        Hdriver += B * np.sin(theta) * sparse_kron_n(operators)
    
    return Hzz, Hdriver


def tilted_field_ising_gs_fidsus(n, theta, B, eps=0.00001):
    """
    Returns T = 0 fid sus for tilted field ising model.
    """
    # Compute the ground state vector for H(B)
    hzz, hdriver = sparse_tilted_ising_field_ham(n, theta, B)
    h = hzz + hdriver
    _, evecs = scipy.sparse.linalg.eigsh(h, k=1, which='SA')
    psi0 = evecs[:, 0]
    
    # Compute the ground state vector for H(B + epsilon)
    eps_h = h + eps * hdriver
    _, evecs_eps = scipy.sparse.linalg.eigsh(eps_h, k=1, which='SA')
    psiPeps = evecs_eps[:, 0]
    
    # Compute dPsi
    dPsi = (psiPeps - psi0) / eps
    
    # Compute chi
    dPsi_dag = dPsi.conj().transpose()
    psi0_dag = psi0.conj().transpose()
    chi = dPsi_dag @ dPsi - dPsi_dag @ psi0 * psi0_dag @ dPsi
    #chi = (np.dot(np.conj(dPsi.T), dPsi) - 
    #       np.dot(np.conj(dPsi.T), psi0) * np.dot(np.conj(psi0.T), dPsi))
    if chi.imag < 1e-8:
        chi = chi.real

    return chi
# %%
# %%
def compute_exact_O_observables(n, beta, tau=0.0, ofname="O.txt"):
    """
    Computes <O> for n qubit system.
    """
    h = load_hamiltonian(n, "/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/H.txt")
    o = load_hamiltonian(n, f"/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/{ofname}")

    exp_h = sp.linalg.expm(-beta * h)
    z = exp_h.trace()
    # basic observables
    o_exp_val = (o @ exp_h).trace() / z # <O>
    osq_exp_val = (o @ o @ exp_h).trace() / z # <O^2>
    
    # <O(\tau)O>
    tau_exp_h = sp.linalg.expm(-tau * h)
    taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
    corr_exp_val = (o @ tau_exp_h @ o @ taubeta_exp_h).trace() / z

    # \int_0^{\beta} <O(\tau)O> \dtau 
    def corr_integrand(tau):
        tau_exp_h = sp.linalg.expm(-tau * h)
        taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
        return (o @ tau_exp_h @ o @ taubeta_exp_h).trace() / z

    Eint_val, Eint_error = scipy.integrate.quad(corr_integrand, 0, beta)
    print(Eint_error)

    # \int_0^{\beta/2} \tau <O(\tau)O> \dtau 
    def tau_corr_integrand(tau):
        tau_exp_h = sp.linalg.expm(-tau * h)
        taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
        return tau * (o @ tau_exp_h @ o @ taubeta_exp_h).trace() / z

    Fint_val, Fint_error = scipy.integrate.quad(tau_corr_integrand, 0, beta / 2)
    print(Fint_error)
    

    return o_exp_val, osq_exp_val, corr_exp_val, Eint_val, Fint_val

def compute_exact_AB_observables(n, beta, tau=0.0, Afname="A.txt", Bfname="B.txt"):
    """
    Computes <O> for n qubit system.
    """
    h = load_hamiltonian(n, "/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/H.txt")
    a = load_hamiltonian(n, f"/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/{Afname}")
    b = load_hamiltonian(n, f"/Users/worknic/research_projects/fidsus_pmrqmc_prl_response_code/private_working_repo/{Bfname}")

    exp_h = sp.linalg.expm(-beta * h)
    z = exp_h.trace()
    # basic observables
    ab_exp_val = (a @ b @ exp_h).trace() / z # <AB>
    
    # <A(\tau)B>
    tau_exp_h = sp.linalg.expm(-tau * h)
    taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
    corr_exp_val = (a @ tau_exp_h @ b @ taubeta_exp_h).trace() / z

    # \int_0^{\beta} <O(\tau)O> \dtau 
    def corr_integrand(tau):
        tau_exp_h = sp.linalg.expm(-tau * h)
        taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
        return (a @ tau_exp_h @ b @ taubeta_exp_h).trace() / z

    Eint_val, Eint_error = scipy.integrate.quad(corr_integrand, 0, beta)
    print(Eint_error)

    # \int_0^{\beta/2} \tau <O(\tau)O> \dtau 
    def tau_corr_integrand(tau):
        tau_exp_h = sp.linalg.expm(-tau * h)
        taubeta_exp_h = sp.linalg.expm((tau-beta) * h)
        return tau * (a @ tau_exp_h @ b @ taubeta_exp_h).trace() / z

    Fint_val, Fint_error = scipy.integrate.quad(tau_corr_integrand, 0, beta / 2)
    print(Fint_error)
    

    return ab_exp_val, corr_exp_val, Eint_val, Fint_val
# %%
