//==============================================================================
// DIVDIFF MODULE - Divided Differences for Exponential Functions
//==============================================================================
//
// OVERVIEW:
// This module implements efficient algorithms for calculating divided differences
// of the exponential function, which are essential for quantum Monte Carlo simulations.
// The implementation uses extended-precision arithmetic to maintain numerical stability
// with extreme exponent values.
//
// KEY ALGORITHMS:
// 1. ExExFloat: Extended-precision floating-point arithmetic for handling extreme exponents
// 2. Divided difference calculation with dynamic scaling for numerical stability
// 3. Efficient addition and removal of input points without full recalculation
//
// MATHEMATICAL BACKGROUND:
// Divided differences of the exponential function are defined as:
//   d[z0,z1,...,zk] = (e^{z1} - e^{z0})/(z1 - z0) for k=1
//   d[z0,z1,...,zk] = (d[z1,...,zk] - d[z0,...,zk-1])/(zk - z0) for k>1
//
// PRIMARY REFERENCES:
// - Gupta, L., Barash, L., Hen, I. (2020). "Calculating the divided differences of the
//   exponential function by addition and removal of inputs." Computer Physics Communications 254, 107385.
// - Ezzell, N., Barash, L., Hen, I. (2024). "Exact and universal quantum Monte Carlo
//   estimators for energy susceptibility and fidelity susceptibility." arXiv:2408.03924.
// - Ezzell, N., Hen, I. (2025). "Advanced measurement techniques in quantum Monte Carlo:
//   The permutation matrix representation approach." arXiv:2504.07295.
//
// USAGE:
// 1. Initialize: divdiff_init() - call once before using any divdiff functions
// 2. Create: divdiff dd(max_points, max_scale) - initialize with capacity
// 3. Add points: dd.AddElement(z_value) - add points incrementally
// 4. Access results: dd.dividedDifferences[i] contains i! * d[z0,...,zi]
// 5. Cleanup: divdiff_clear_up() - call once when done (optional)
//
// PERFORMANCE:
// - Addition: O(n) per point where n is current number of points
// - Removal: O(n) per point
// - Memory: O(n * max_scale) for storing intermediate results
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//==============================================================================

#pragma once

#include<random>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<vector>

//==============================================================================
// GLOBAL CONSTANTS AND VARIABLES
//==============================================================================

// Lookup table for inverse powers of 2 (1/2^i) for extended-precision arithmetic
// Pre-computed for efficiency in ExExFloat operations
extern double* invPow2;

// Maximum exponent value supported by the extended-precision arithmetic
// Determines the range of exponents that can be handled without overflow
extern const int MAX_EXP_RANGE;

// Extra length allocated for arrays to prevent boundary issues
// Provides buffer space for Hermite polynomial calculations
extern const int EXTRA_LEN;

//==============================================================================
// UTILITY FUNCTIONS
//==============================================================================

// Template functions for finding minimum and maximum values
template <typename T> T min (T a, T b) { return (b>=a)?a:b;}
template <typename T> T max (T a, T b) { return (b>=a)?b:a;}

//==============================================================================
// EXEXTENDED-PRECISION FLOATING POINT CLASS
//==============================================================================

/**
 * Extended-precision floating-point arithmetic class for handling extreme exponent values.
 *
 * This class represents numbers in the form: mantissa * 2^exponent
 * where mantissa is in [0.5, 1.0) and exponent is an integer.
 * This allows handling very large and very small numbers that would overflow/underflow
 * standard double precision, which is essential for quantum Monte Carlo calculations.
 *
 * Key features:
 * - Automatic normalization after each operation
 * - Support for all basic arithmetic operations
 * - Efficient comparison and conversion functions
 * - Square root operation with proper exponent handling
 */
class ExExFloat{
private:
    double mantissa;      // Normalized mantissa in [0.5, 1.0)
    int exponent;         // Integer exponent for 2^exponent scaling

public:
    // Constructors
    ExExFloat();                                              // Default constructor (0.5 * 2^1 = 1.0)
    ExExFloat(double value);                                  // Construct from double (auto-normalizes)
    ExExFloat(const ExExFloat& other);                        // Copy constructor

    // Core operations
    void normalize();                                          // Normalize mantissa to [0.5, 1.0)
    void initExpMu(double mu);                                // Initialize as e^mu for scaling
    void print() const;                                        // Print value (scientific notation for large exponents)

    // Assignment operators
    ExExFloat operator =(const ExExFloat& other);             // Assignment from ExExFloat
    ExExFloat operator =(double value);                        // Assignment from double

    // Arithmetic operators
    ExExFloat operator +(const ExExFloat& other) const;       // Addition
    ExExFloat operator -(const ExExFloat& other) const;       // Subtraction
    ExExFloat operator +=(const ExExFloat& other);            // Addition assignment
    ExExFloat operator -=(const ExExFloat& other);            // Subtraction assignment
    ExExFloat operator *(const ExExFloat& other) const;       // Multiplication
    ExExFloat operator /(const ExExFloat& other) const;       // Division
    ExExFloat operator *=(const ExExFloat& other);            // Multiplication assignment
    ExExFloat operator /=(const ExExFloat& other);            // Division assignment

    // Scalar operations
    ExExFloat operator *(double value) const;                 // Multiplication by scalar
    ExExFloat operator /(double value) const;                 // Division by scalar
    ExExFloat operator *=(double value);                       // Scalar multiplication assignment
    ExExFloat operator /=(double value);                       // Scalar division assignment

    // Friend functions for reverse operations
    friend ExExFloat operator *(double lhs, const ExExFloat& rhs);  // Double * ExExFloat
    friend ExExFloat operator /(double lhs, const ExExFloat& rhs);  // Double / ExExFloat

    // Comparison operators
    int operator >=(double value) const;                       // Greater than or equal to double
    int operator >=(const ExExFloat& other) const;            // Greater than or equal to ExExFloat

    // Utility functions
    double getDouble() const;                                  // Convert to standard double
    int sign() const;                                          // Get sign (-1, 0, +1)
    ExExFloat absoluteValue() const;                          // Absolute value
    ExExFloat squareRoot() const;                             // Square root with proper exponent handling
};

//==============================================================================
// GLOBAL INITIALIZATION FUNCTIONS
//==============================================================================

/**
 * Initialize the global lookup table for extended-precision arithmetic.
 *
 * This function must be called once before using any divdiff functionality.
 * It pre-computes the inverse powers of 2 lookup table for efficient ExExFloat operations.
 *
 * @note Call this function once at program startup.
 */
void divdiff_init();

/**
 * Clean up global resources used by the divdiff module.
 *
 * This function frees the memory allocated for the inverse powers of 2 lookup table.
 *
 * @note Call this function once at program cleanup (optional but recommended).
 */
void divdiff_clear_up() noexcept;

//==============================================================================
// DIVIDED DIFFERENCES CALCULATOR CLASS
//==============================================================================

/**
 * Main class for calculating divided differences of the exponential function.
 *
 * This class implements efficient algorithms for computing divided differences
 * d[z0,z1,...,zk] of the exponential function exp(z), which are fundamental
 * for quantum Monte Carlo simulations. The implementation supports dynamic
 * addition and removal of points without full recalculation, and uses
 * extended-precision arithmetic for numerical stability.
 *
 * Key algorithms:
 * - Hermite polynomial expansion around a central point
 * - Dynamic scaling for numerical stability
 * - Efficient O(n) addition and removal operations
 *
 * Usage:
 *   divdiff dd(1000, 100);  // Support up to 1000 points, scaling factor 100
 *   dd.AddElement(1.0);     // Add first point
 *   dd.AddElement(2.0);     // Add second point
 *   // dd.dividedDifferences[i] contains i! * d[z0,...,zi]
 */
class divdiff{

protected:
    // Internal arrays for computation
    double* zBackups;                                // Backup storage for z array during reallocation
    ExExFloat* hCoeffs;                             // Hermite polynomial coefficients h_k
    ExExFloat* derivTerms;                          // Storage for derivative terms d[k][n]

    // Algorithm parameters
    int scaleFactor;                                 // Current scaling factor for numerical stability
    int maxLen;                                    // Maximum length of input arrays (default: 10001)
    int maxScale;                                   // Maximum scaling factor supported (default: 500)

    // Central expansion point
    double meanVal;                             // Central value μ for Hermite expansion
    ExExFloat expMu;                            // e^μ for efficient calculations

    // Helper functions
    double calculateMean(double* values, int count);                // Calculate arithmetic mean
    double calculateMaxAbsDiff(double* values, int length); // Find max |z_i - z_j|

public:
    // Public data arrays
    double* zVals;                                  // Array of input z values
    ExExFloat* divDiffs;                            // Results: i! * d[z0,...,zi]
    int currentLength;                              // Current number of points stored

    // Constructors and destructor
    divdiff(int maxLen_, int maxScale_);                // Constructor
    divdiff(const divdiff& other);                          // Copy constructor
    divdiff& operator=(const divdiff& other);               // Copy assignment
    ~divdiff();                                             // Destructor

    // Memory management
    void allocateMemory();                                // Allocate internal arrays
    void freeMemory();                                    // Free internal arrays
    long long int getMemoryUsage();                       // Get total memory usage in bytes

    // Debug and utility functions
    void printExExFloatList(ExExFloat* list, int length, const char* name);
    void printDoubleList(double* list, int length, const char* name);
    void printIntegerList(int* list, int length, const char* name);
    void printVectorList(std::vector<int> list, const char* name);

    // Array backup/restore for reallocation
    void backupZValues(int length);                      // Backup z array to temporary storage
    void restoreZValues(int length);                     // Restore z array from backup

    // Scaling factor management
    int scalingChanged();                                  // Check if scaling factor needs adjustment
    void adjustScalingFactorIfNeeded();                  // Adjust scaling if needed

    // Core divided difference operations
    void addElement(double newZValue, int forcedScale = 0, double forcedCentral = 0);
    void removeElement();                               // Remove most recently added element
    int removeValue(double targetValue, int forcedScale = 0, double forcedCentral = 0);
    void addAllElements(int length, int forcedScale = 0);  // Bulk addition of elements
};

using DivDiff = divdiff;

