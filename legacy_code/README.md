-----------------------------------------------------------------------------------------------------------

This repository contain code to run simulations contained in the preprint: Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility (2024). 

Link to asssociated code Zenodo (like this repo but fixed in time as when published, i.e., no commits): <a href="https://zenodo.org/records/13251139">https://zenodo.org/records/13251139</a>

Link to associated data/plotting Zenodo: <a href="https://zenodo.org/records/13251202">https://zenodo.org/records/13251202</a>

-----------------------------------------------------------------------------------------------------------
# Attribution
Most of the simulation code was (privately) forked from <a href="https://github.com/LevBarash/PMRQMC">PMRQMC</a>, as developed in: Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians, Physical Review Research 6, 013281 (2024).

# Instructions
Details on using PMR-QMC itself are contained in the aforementioned <a href="https://github.com/LevBarash/PMRQMC">PMRQMC</a>. To replicate the experiments in [preprint], simply run the python driver files contained in the "x_experiments/" folders. At minimum, you will need a version of python that contains `numpy`, but one containing also `scipy` allows one to generate exact results for comparison as well. The fidsus_driver takes 3 possible command line arguments... 
- `python fidsus_experiments/zhang_prl_model_fidsus_driver.py 0` runs a simple test to ensure things are working properly for susceptibility simulations
- `python fidsus_experiments/zhang_prl_model_fidsus_driver.py 1` runs simulation to generate data for Figure 1
- `python fidsus_experiments/zhang_prl_model_fidsus_driver.py 2` runs simulation to generate data for Figure 2c and 2d

In each case, several temporary simulation files will be created in the base directory, and an associated data file (as a csv) will be placed in the `fidsus_experiments/` directory.

Similarly, 
- `python correlator_experiments/zhang_prl_model_correlator_driver.py 0` runs a simple test to ensure things are working properly for correlator simulations
- `python correlator_experiments/zhang_prl_model_correlator_driver.py 2` runs simulation to generate data for Figure 2a and 2b

# Brief comments on code structure and function
- Code in the base directory is forked from <a href="https://github.com/LevBarash/PMRQMC">PMRQMC</a>. The code in `mainqmc.hpp` computes the susceptibility estimators and that in `correlator_mainqmc.hpp` computes the correlator estimators. (Technically, they both can compute either, but for technical reasons, it was easier to make separate files for each).
- `Utils/` contains new code written entirely for [preprint]. Each file has a brief description of its purpose at the top.
- Both `fidsus_experiments/` and `correlator_experiments/` just contain simulation drivers as explained above
