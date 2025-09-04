
//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians, Physical Review Research 6, 013281 (2024).
//
// Various advanced measurement capabilities were added as part of the
// work introduced in the papers:
// * Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
// * Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).
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
    
#define Tsteps 1000 // number of Monte-Carlo initial equilibration updates
#define steps 10000 // number of Monte-Carlo updates
#define stepsPerMeasurement 1 // number of Monte-Carlo updates per measurement
#define beta 0.1 // inverse temperature
#define tau 0.05 //imaginary propogation time
#define parity_cond 0 // controls parity subspace measurement 

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
#define MEASURE_HDIAG_EINT
#define MEASURE_HDIAG_FINT
#define MEASURE_HOFFDIAG_CORR
#define MEASURE_HOFFDIAG_EINT
#define MEASURE_HOFFDIAG_FINT

//
// Below are the implementation parameters:
//

#define qmax     1000                // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
#define EXHAUSTIVE_CYCLE_SEARCH      // comment this line for a more restrictive cycle search
#define GAPS_GEOMETRIC_PARAMETER 0.8 // parameter of geometric distribution for the length of gaps in the cycle completion update
#define COMPOSITE_UPDATE_BREAK_PROBABILITY  0.9   // exit composite update at each step with this probability

// #define ABS_WEIGHTS                  // uncomment this line to employ absolute values of weights rather than real parts of weights
// #define EXACTLY_REPRODUCIBLE         // uncomment this to always employ the same RNG seeds and reproduce exactly the same results

//
// Uncomment or comment the macros below to enable or disable the ability to checkpoint and restart
//
