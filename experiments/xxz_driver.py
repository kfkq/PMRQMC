"""
Suscepibility experiment driver
for XXZ model
H = (Jz/4) \sum_{<i,j>}Z_iZ_j + (lam/4) \sum_{<i,j>} (X_iX_j + Y_iY_j)
and random rotations thereof. The division by four is from the fact
original model uses Sz = (1/2) Z, but we use direct Pauli matrices.

Supports 1D chain, 2D square lattice,
or 2D triangular lattice with or
without periodic boundary conditions (PBC)

This driver was written in support of the experiments carried out in:
// * Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
// * Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).
"""
# %%
import subprocess, datetime, sys, os
import numpy as np
sys.path.append("..")
sys.path.append("../utils")
from pauli_manipulations import PauliTerm, PauliH, PauliU
from ioscripts import make_all_stand_param_fstr
from build_pauliH_recipes import (build_1d_xxz_xxyy, build_1d_xxz_zz,
                                  build_square_xxz_xxyy, build_square_xxz_zz,
                                  build_triangle_xxz_xxyy, build_triangle_xxz_zz)

def main(nt, L, lat, pbc, lam, beta, tau, parity, strnow, eps=None, l=None, Tsteps=1000000, steps=10000000, stepsPerMeasurement=10,
         save=True, restart=False, diag_pert=True):
    # ==============================================
    # Choose random parametrs if not given 
    # ==============================================
    if lat == "chain":
        n = L
    else:
        n = L * L
    if l is None:
        l = np.random.choice(range(1, 2*n + 2))
        
    # make relevant directory name, copy over important files
    #dir_name = f"../{lat}_tfim/square_tfim_L_{L}_lam_{lam}"
    dir_name = ".."
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    # ==============================================
    # Building H.txt
    # ==============================================
    if diag_pert is True:
        if lat == "chain":
            h0 = build_1d_xxz_xxyy(L, pbc)
            h1 = build_1d_xxz_zz(L, pbc)
        elif lat == "square":
            h0 = build_square_xxz_xxyy(L, pbc)
            h1 = build_square_xxz_zz(L, pbc)
        elif lat == "triangle":
            h0 = build_triangle_xxz_xxyy(L, pbc)
            h1 = build_triangle_xxz_zz(L, pbc)
    # just flip h1/h0 if we want off-diag driver
    else:
        if lat == "chain":
            h1 = build_1d_xxz_xxyy(L, pbc)
            h0 = build_1d_xxz_zz(L, pbc)
        elif lat == "square":
            h1 = build_square_xxz_xxyy(L, pbc)
            h0 = build_square_xxz_zz(L, pbc)
        elif lat == "triangle":
            h1 = build_triangle_xxz_xxyy(L, pbc)
            h0 = build_triangle_xxz_zz(L, pbc)

    # rotate the Hamiltonian
    u = PauliU(n)
    if l > 0:
        u.set_as_random(l, eps)
        uh0 = h0.conjugate(u)
        uh1 = h1.conjugate(u)
    else:
        uh0 = h0
        uh1 = h1
    uh = uh0 + lam * uh1
    # save Hamiltonian file
    with open(f"{dir_name}/H.txt", 'w') as f:
        f.write(uh.to_pmr_str())

    # ==============================================
    # Building observables A and B
    # ==============================================
    # A = Z1Z2
    pzz = PauliTerm(1.0, [1, 2], ['Z','Z'], n)
    a = PauliH(n, [pzz])
    # B = X1X2 + Y1Y2
    pxx = PauliTerm(1.0, [1, 2], ['X','X'], n)
    pyy = PauliTerm(1.0, [1, 2], ['Y', 'Y'], n)
    b = PauliH(n, [pxx, pyy])
    with open(f"{dir_name}/A.txt", 'w') as f:
        f.write(a.to_pmr_str())
    with open(f"{dir_name}/B.txt", 'w') as f:
        f.write(b.to_pmr_str())
    # ==============================================
    # Running PMR QMC
    # ==============================================
    param_str = make_all_stand_param_fstr(beta, tau, Tsteps, steps, stepsPerMeasurement, parity, save, restart)
    with open(f"{dir_name}/parameters.hpp", "w") as f:
        f.write(param_str)
    # compile and run PMR-QMC
    job_file = "job_file.sh"
    temp_out_file = f"{lat}_xxz_experiment_{strnow}.txt"
    with open(f"{dir_name}/{job_file}", 'w') as fh:
        fh.write("#!/bin/bash\n")
        # compile and run
        fh.write("g++ -O3 -std=c++11 -o prepare.bin prepare.cpp\n")
        fh.write("./prepare.bin H.txt A.txt A.txt A.txt A.txt A.txt B.txt B.txt B.txt B.txt B.txt A.txt A.txt A.txt A.txt A.txt\n")
        if nt > 1:
            fh.write("mpicxx -O3 -std=c++11 -o PMRQMC_mpi.bin PMRQMC_mpi.cpp\n")
            #fh.write(f"mpirun -n {nt} ./PMRQMC_mpi.bin > {temp_out_file}\n")
            fh.write(f"mpirun -n {nt} ./PMRQMC_mpi.bin")
        else:
            fh.write("g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp\n")
            #fh.write(f"./PMRQMC.bin > {temp_out_file}")
            fh.write(f"./PMRQMC.bin")
    subprocess.run(f"cd {dir_name}; chmod +x {job_file}; ./{job_file}", shell=True)
    print(f"Successfully submitted {job_file}")

    return
# %%

# %%
if __name__=="__main__":
    # get current date-time
    now = datetime.datetime.now()
    strnow = now.strftime("%Y-%m-%d_%H-%M-%S")
    # hard-coding input parameters
    # ======================================
    # basic model parameters
    # ======================================
    # lat -- lattice type: chain, square, triangle
    lat = "square"
    # L -- lattice size: L spins for chain, LxL for square, triangle
    L = 2
    # pbc -- periodic boundary conditions: 0 is False 1 is True
    pbc = 0
    # lam -- transverse field strength: a float
    lam = 0.1
    # ======================================
    # Basic simulation parameters
    # ======================================
    # beta -- inverse temperature: positive float
    beta = 0.1
    # Tsteps -- equilibration steps: positive int
    Tsteps = int(1e3)
    # steps -- QMC updates: positive int
    steps = int(1e4)
    # stepsPer... -- : # QMC updates per measurement: positive int
    stepsPerMeasurement = int(1)
    # ======================================
    # Simulation meta-parameters
    # ======================================
    # save -- if True, saves QMC simulation state at end of calculation: Bool
    save = False
    # restart -- if True, restarts calculation where left off: Bool
    restart = False
    # nt -- number of threads to run with MPI: positive int
    # for nt > 1, may fail due to OSX permissions. If so,
    # simply execute ./job_file.sh directly in base directory.
    nt = 1
    # ======================================
    # Advanced parameters 
    # ======================================
    # tau -- imaginary time: float in range [0.0, beta]
    tau = beta / 2
    # diag_pert -- if True, driving term is diagonal (ZZ term)
    diag_pert = False
    # l -- number of terms in random Pauli unitary: non-negative int
    # rotates Hamiltonian if l > 0 (see Appendix A, arXiv:2408.03924v1)
    l = 0 
    # ======================================
    # Set-up and execute simulation
    # ======================================
    # the hard-coded 0 is for parity, which is only a relevant setting for the TFIM for now (see Appendix A, arXiv:2408.03924v1)
    main(nt, L, lat, pbc, lam, beta, tau, 0, strnow, l=l, Tsteps=Tsteps, steps=steps, stepsPerMeasurement=stepsPerMeasurement,save=save, restart=restart, diag_pert=diag_pert)
 # %%
