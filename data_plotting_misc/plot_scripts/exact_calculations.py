"""
//
// This program contains code to do "exact numerical calculations" in support of:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
//
//
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
    result, error = scipy.integrate.quad(lambda tau: 2 * prl_gs_gtau(Bz, tau), 0, np.inf)
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
