This directory contains input driver files in relation to the preprints:
[1] Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
[2] Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).

-----------------------------------------------------------------------------------------------------------

Data in the subdirectory `tfim_observables' are in relation to [2] ([arXiv:2504.07295]). These contain text file descriptions that can be fed to our QMC algorithm for the custom observables defined in Eqs. (10)--(13) of the manuscript.

Data in `xxz_H_files' are in relation to [1] ([arXiv:2408.03924]). This contains the H.txt files that define the XXZ model on a square lattice with periodic boundary conditions of various sizes. This is not necessary (one can use our XXZ driver), but it speeds up pre-processing to not have to regenerate these on the fly.