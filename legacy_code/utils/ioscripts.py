"""
//
// This program contains code to handle common I/O tasks in support of:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
//
//
"""
# %%
import re
import numpy as np

def make_no_stand_param_fstr(seed, beta, tau, Tsteps=100000, steps=1000000,
stepsPerMeasurement=10):
    """
    Writes parameters.hpp without any standard observables.
    """
    param_str = make_param_file_header(seed, beta, tau, Tsteps, steps, stepsPerMeasurement)
    param_str += make_param_file_footer()
    return param_str

def make_hdiag_fidsus_param_fstr(seed, beta, Tsteps=100000, steps=1000000,
stepsPerMeasurement=10):
    """
    Makes parameters.hpp with Hdiag measurements relevant to compute
    fidelity suscepibility.
    """
    param_str = make_param_file_header(seed, beta, 0.0, Tsteps, steps, stepsPerMeasurement)
    param_str += """
//
// Below is the list of standard observables:
//

#define MEASURE_HDIAG                // <H_{diag}>      is measured when this line is not commented
#define MEASURE_HDIAG_INT1
#define MEASURE_HDIAG_INT2
#define MEASURE_HDIAG_INT3
"""
    param_str += make_param_file_footer()

    return param_str

def make_hoffdiag_fidsus_param_fstr(seed, beta, Tsteps=100000, steps=1000000,
stepsPerMeasurement=10):
    """
    Makes parameters.hpp with Hoffdiag measurements relevant to compute
    fidelity suscepibility.
    """
    param_str = make_param_file_header(seed, beta, 0.0, Tsteps, steps, stepsPerMeasurement)
    param_str += """
//
// Below is the list of standard observables:
//

#define MEASURE_HOFFDIAG                // <H_{diag}>      is measured when this line is not commented
#define MEASURE_HOFFDIAG_INT1
#define MEASURE_HOFFDIAG_INT2
#define MEASURE_HOFFDIAG_INT3
"""
    param_str += make_param_file_footer()

    return param_str

def make_all_stand_param_fstr(seed, beta, tau, Tsteps=100000, steps=1000000,
stepsPerMeasurement=10):
    """
    Makes parameters.hpp with all standard observables.
    """
    param_str = make_param_file_header(seed, beta, tau, Tsteps, steps, stepsPerMeasurement)
    param_str += """
//
// Below is the list of standard observables:
//

#define MEASURE_H                    // <H>             is measured when this line is not commented
#define MEASURE_H2                   // <H^2>           is measured when this line is not commented
#define MEASURE_HDIAG                // <H_{diag}>      is measured when this line is not commented
#define MEASURE_HDIAG2               // <H_{diag}^2>    is measured when this line is not commented
#define MEASURE_HOFFDIAG             // <H_{offdiag}>   is measured when this line is not commented
#define MEASURE_HOFFDIAG2            // <H_{offdiag}^2> is measured when this line is not commented
#define MEASURE_Z_MAGNETIZATION      // Z-magnetization is measured when this line is not commented
#define MEASURE_HDIAG_CORR
#define MEASURE_HDIAG_INT1
#define MEASURE_HDIAG_INT2
#define MEASURE_HDIAG_INT3
#define MEASURE_HOFFDIAG_CORR
#define MEASURE_HOFFDIAG_INT1
#define MEASURE_HOFFDIAG_INT2
#define MEASURE_HOFFDIAG_INT3
"""
    param_str += make_param_file_footer()

    return param_str

def parse_otxt_temp(temp_fname):
    """
    Parse temporary output file.
    """
    with open(temp_fname, 'r') as f:
        lines = f.readlines()
    # define regular expression parser for numbers
    #p = '[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+'
    p = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|nan'
    # get emergent quantities
    sign = re.findall(p, lines[2])[0]
    sign_std = re.findall(p, lines[3])[0]
    mean_q = re.findall(p, lines[4])[0]
    max_q = re.findall(p, lines[5])[0]
    emergent = [sign, sign_std, mean_q, max_q]
    otxt_array = []
    otxt_std = []
    # get O.txt observables
    for j in range(1, 7):
        idx = 6 + 3*j - 2
        otxt_array.append(re.findall(p, lines[idx])[0])
        idx += 1
        otxt_std.append(re.findall(p, lines[idx])[0])
    otxt_array = np.array([float(x) for x in otxt_array])
    otxt_std = np.array([float(x) for x in otxt_std])
    # get time simulation took
    idx = 6 + 3*6
    time = re.findall(p, lines[idx])[0]

    return emergent, otxt_array, otxt_std, time

def parse_fidsus_temp(temp_fname):
    """
    Parse temporary output file.
    """
    with open(temp_fname, 'r') as f:
        lines = f.readlines()
    # define regular expression parser for numbers
    #p = '[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+'
    p = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|nan'
    # get emergent quantities
    sign = re.findall(p, lines[2])[0]
    sign_std = re.findall(p, lines[3])[0]
    mean_q = re.findall(p, lines[4])[0]
    max_q = re.findall(p, lines[5])[0]
    emergent = [sign, sign_std, mean_q, max_q]
    obs_array = []
    obs_std = []
    # get <O> (or Hdiag, Hoffdiag), and 3 integrals thereof
    for j in range(1, 5):
        idx = 6 + 3*j - 2
        obs_array.append(re.findall(p, lines[idx])[0])
        idx += 1
        obs_std.append(re.findall(p, lines[idx])[0])
    # get time simulation took
    idx = 6 + 3*4
    time = re.findall(p, lines[idx])[0]

    return emergent, obs_array, obs_std, time

def parse_correlator_temp(temp_fname):
    """
    Parse temporary output file.
    """
    with open(temp_fname, 'r') as f:
        lines = f.readlines()
    # define regular expression parser for numbers
    #p = '[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+'
    p = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|nan'
    # get emergent quantities
    sign = re.findall(p, lines[2])[0]
    sign_std = re.findall(p, lines[3])[0]
    mean_q = re.findall(p, lines[4])[0]
    max_q = re.findall(p, lines[5])[0]
    emergent = [sign, sign_std, mean_q, max_q]
    obs_array = []
    obs_std = []
    # get <O> and <O(\tau)O>
    for j in range(1, 3):
        idx = 6 + 3*j - 2
        obs_array.append(re.findall(p, lines[idx])[0])
        idx += 1
        obs_std.append(re.findall(p, lines[idx])[0])
    # get time simulation took
    idx = 6 + 3*2
    time = re.findall(p, lines[idx])[0]

    return emergent, obs_array, obs_std, time

#
# helper functions only for this file
#
def make_param_file_header(seed, beta, tau, Tsteps=10000000, steps=1000000,
stepsPerMeasurement=10):
    """
    Writes parameters.hpp header.
    """
    param_str = """
//
// This program helps implement estimators for the titular quantities in the preprint:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
// 
// This work augments the base code whose information is below.
//
//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians, Physical Review Research 6, 013281 (2024).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

//
// Below are the parameter values:
//
    """
    # add actual adjustable parameters
    param_str += f"""
#define Tsteps {Tsteps} // number of Monte-Carlo initial equilibration updates
#define steps {steps} // number of Monte-Carlo updates
#define stepsPerMeasurement {stepsPerMeasurement} // number of Monte-Carlo updates per measurement
#define beta {beta} // inverse temperature
#define tau {tau} //imaginary propogation time
#define rng_seed {seed}
"""
    return param_str

def make_param_file_footer():
    # add final technical parameters, not usually needed to adjust
    param_str = """
//
// Below are the implementation parameters:
//

#define qmax     1000                // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
#define EXHAUSTIVE_CYCLE_SEARCH      // comment this line for a more restrictive cycle search
    """
    return param_str
#%%
