"""
//
// This auxillary program builds, runs, and collects data from a QMC estimator code for the preprint:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
//
//
Driver code to estimate energy susceptibility and fidelity susceptibility
for 2 qubit PRL model in [10.1103/PhysRevLett.100.100501] and special
rotations of this model described in [preprint] mentioned above.
"""
# %%
import subprocess, datetime, sys, os
import numpy as np
sys.path.append("..")
sys.path.append("../utils")
from pauli_manipulations import PauliTerm, PauliH, PauliU
from ioscripts import (make_no_stand_param_fstr, parse_fidsus_temp)

def main(n, seed, lam, beta, strnow, eps=None, l=None, Tsteps=100000, steps=1000000, stepsPerMeasurement=10):
    # ==============================================
    # Choose random parametrs if not given 
    # ==============================================
    np.random.seed(seed)
    if l is None:
        l = np.random.choice(range(1, 2*n + 2))
    # ==============================================
    # Building H.txt
    # ==============================================
    if n < 2:
        raise ValueError("n must be greater than or equal to 2")
    # define basic 2q model from PRL Hamiltonian
    x1 = PauliTerm(1.0, [1], ['X'], n)
    x2 = PauliTerm(1.0, [2], ['X'], n)
    z1 = PauliTerm(1.0, [1], ['Z'], n)
    z2 = PauliTerm(1.0, [2], ['Z'], n)
    z1z2 = PauliTerm(1.0, [1, 2], ['Z', 'Z'], n)
    h0 = PauliH(n, [z1z2, 0.1 * x1, 0.1 * x2])
    h1 = PauliH(n, [z1, z2])

    # rotate the Hamiltonian
    u = PauliU(n)
    if l > 0:
        u.set_as_random(l, eps, seed=seed)
    uh0 = h0.conjugate(u)
    uh1 = h1.conjugate(u)
    uh = uh0 + lam * uh1

    # ==============================================
    # Running PMR QMC
    # ==============================================
    # save Hamiltonian file and O.txt file
    with open("../H.txt", 'w') as f:
        f.write(uh.to_pmr_str())
    with open("../O.txt", 'w') as f:
        f.write(uh1.to_pmr_str())
    param_str = make_no_stand_param_fstr(seed, beta, 0.0, Tsteps, steps, stepsPerMeasurement)
    with open("../parameters.hpp", "w") as f:
        f.write(param_str)
    # compile and run PMR-QMC
    job_file = "job_file.sh"
    temp_out_file = "temp_fidsus_data_zhang_prl_rot_model.txt"
    with open("../" + job_file, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("g++ -O3 -std=c++11 -o prepare.bin prepare.cpp\n")
        fh.write("./prepare.bin H.txt $(ls O.txt O.txt O.txt O.txt  2> /dev/null)\n")
        #fh.write("./prepare.bin H.txt\n")
        fh.write("g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp\n")
        fh.write(f"./PMRQMC.bin > {temp_out_file}")
    print(f"Trying to submit: {job_file}")
    subprocess.run(f"cd ..; chmod +x {job_file}; ./{job_file}", shell=True)
    #subprocess.run(["bash", job_file])
    print("Success?")
    # ==============================================
    # Process and save output
    # ==============================================
    emergent, obs, std, time = parse_fidsus_temp("../" + temp_out_file)
    #----- form file ------
    fname = f"fidsus_data_zhang_prl_rot_curve_{strnow}.csv"
    # if file does not exit, open and add header
    if not os.path.exists(fname):
        with open(fname, "w") as f:
            header1 = "n,lam,beta,H1,H1_std,"
            header2 = "int1,int1_std,"
            header3 = "int2,int2_std,"
            header4 = "int3,int3_std,"
            header5 = "rng,eps,l,Tsteps,steps,stepsPerMeasurement,"
            header6 = "sign,sign_std,q,qmax,time"
            f.write(header1 + header2 + header3 + header4 + header5 + header6)
    # always append
    with open(fname, "a+") as f:
        # save plot param inputs
        line1 = f"\n{n},{lam},{beta},{obs[0]},{std[0]},"
        line2 = f"{obs[1]},{std[1]},"
        line3 = f"{obs[2]},{std[2]},"
        line4 = f"{obs[3]},{std[3]},"
        line5 = f"{seed},{eps},{l},{Tsteps},{steps},{stepsPerMeasurement},"
        line6 = f"{emergent[0]},{emergent[1]},{emergent[2]},{emergent[3]},{time}"
        f.write(line1 + line2 + line3 + line4 + line5 + line6)

    return fname
# %%

# %%
if __name__=="__main__":
    # ====================================
    # Header part of main 
    # ====================================
    # arg parse which experiment to do (see below)
    experiment = int(sys.argv[1])
    # get current date-time
    now = datetime.datetime.now()
    strnow = now.strftime("%Y-%m-%d_%H-%M-%S")
    # get rng seeds
    with open("../utils/rng_seeds.txt", "r") as f:
        seeds = f.readline().split(',')
    seeds = [int(s) for s in seeds]
    # ====================================
    # Experiments one can run
    # ====================================
    # *********************
    # example experiment
    # *********************
    if experiment == 0:
        n = 2
        rng = seeds[0]
        lam = 0.1
        beta = 0.25
        main(n, rng, lam, beta, strnow)
    # *********************
    # Fig. 1 experiment
    # *********************
    elif experiment == 1:
        n = 2
        for rng in seeds[0:5]:
            for beta in [5.0, 10.0, 15.0, 20.0]:
                for lam in np.linspace(-1.5, 1.5, 200):
                    main(n, rng, lam, beta, strnow, l = 0)
    # *********************
    # Fig. 2 (c),(d) experiments
    # *********************
    elif experiment == 2:
        # get rng seeds that had high sign
        n = 100
        beta = 20.0
        high_sign_seeds = [3026438146, 2781301355, 4206712692, 2270171661, 2834120170, 1370102282, 2172515818,  185933287, 4109413259, 3629902865]
        for rng in high_sign_seeds:
            for lam in np.linspace(-1.5, 1.5, 200):
                main(n, rng, lam, beta, strnow)
 # %%
