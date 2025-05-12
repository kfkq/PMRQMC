This directory contains data wrangling scripts in relation to the preprints:
[1] Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
[2] Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).

----------------------------------------------------------------------------------------------------------

Scripts with `advmeas' are relevant to [2] ([arXiv:2504.07295]) and otherwise in relation to [1] ([arXiv:2408.03924]).

These scripts were used to wrangle and compile data we obtained from the Discovery USC HPC cluster. Using a custom Python script to prepare and submit jobs (not shared, as it is specific to our cluster), the output of our QMC simulations were output to text files in many subdirectories with parameter metadata as their names. By combing through these, we extracted the useful simulation parameters and outputs using the `group_..._data.sh` files. We then combined all the outputs into a csv file using the `form_..._csv_summary.sh` script. We include these to give some idea about how one might ramp up simple single QMC scripts with text file outputs to large scale simulations with many different inputs.
