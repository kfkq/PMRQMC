//==============================================================================
// DIVDIFF IMPLEMENTATION - Divided Differences for Exponential Functions
//==============================================================================
//
// Implementation of efficient algorithms for calculating divided differences
// of the exponential function with extended-precision arithmetic support.
//
// References:
// - Gupta, L., Barash, L., Hen, I. (2020). Computer Physics Communications 254, 107385.
// - Ezzell, N., Barash, L., Hen, I. (2024). arXiv:2408.03924.
// - Ezzell, N., Hen, I. (2025). arXiv:2504.07295.
//
// License: Creative Commons Attribution 4.0 International License
//==============================================================================

#include <divdiff.h>
#include <cmath>

//==============================================================================
// GLOBAL CONSTANTS AND VARIABLES
//==============================================================================

// Global lookup table for inverse powers of 2 (1/2^i)
// Pre-computed for efficient extended-precision arithmetic operations
double* invPow2 {NULL};

// Maximum exponent range supported by extended-precision arithmetic
// Limits the range of exponents that can be handled without overflow
const int MAX_EXP_RANGE = 100000;

// Extra array length for boundary safety
// Provides buffer space for Hermite polynomial coefficient calculations
const int EXTRA_LEN = 10;

//==============================================================================
// NUMERICAL CONSTANTS
//==============================================================================

// Conversion factors for logarithm base changes
const double LOG2_E = 1.4426950408889634;    // log2(e) for exponent calculations
const double LOG10_2 = 0.3010299956639812;  // log10(2) for decimal output formatting

// Scaling threshold for numerical stability
const double SCAL_THRESHOLD = 3.5;       // Threshold for scaling factor adjustment

// Default constructor: creates 1.0 (0.5 * 2^1)
ExExFloat::ExExFloat() : mantissa(0.5), exponent(1) {}
// Construct from double: auto-normalize to mantissa * 2^exponent form
ExExFloat::ExExFloat(double obj) : mantissa(obj), exponent(0) { normalize(); }
// Copy constructor: duplicate mantissa and exponent
ExExFloat::ExExFloat(ExExFloat const &obj) : mantissa(obj.mantissa), exponent(obj.exponent) {}
//==============================================================================
// EXEXTENDED-PRECISION FLOAT IMPLEMENTATION
//==============================================================================

// Normalize: adjust mantissa to [0.5, 1.0) and update exponent
void ExExFloat::normalize(){
    int tempExponent;
    mantissa = frexp(mantissa, &tempExponent);
    exponent += tempExponent;
}

// Initialize as e^mu for efficient exponential calculations
void ExExFloat::initExpMu(double mu){
    double exponentValue = mu * LOG2_E;
    exponent = ceil(exponentValue);
    mantissa = pow(2.0, exponentValue - ceil(exponentValue));
}

// Print value in scientific notation for large exponents
void ExExFloat::print() const {
    double exponentBase10, mantissaValue;

    // Convert to base 10 exponent for display
    exponentBase10 = exponent * LOG10_2;
    mantissaValue = mantissa * pow(10, exponentBase10 - floor(exponentBase10));
    exponentBase10 = floor(exponentBase10);

    // Ensure mantissa is at least 1.0 for proper scientific notation
    if(fabs(mantissaValue) < 1.0){
        exponentBase10--;
        mantissaValue *= 10.0;
    }

    // Use standard notation for reasonable exponents, scientific otherwise
    if((exponentBase10 < 99) && (exponentBase10 > -99)){
        printf("%.17f", getDouble());
    } else {
        printf("%.17fe%.0f", mantissaValue, exponentBase10);
    }
}
//==============================================================================
// ASSIGNMENT OPERATORS
//==============================================================================

// Assignment operator: copy mantissa and exponent
ExExFloat ExExFloat::operator =(const ExExFloat& other){
    mantissa = other.mantissa;
    exponent = other.exponent;
    return *this;
}

// Assignment from double: normalize to mantissa * 2^exponent
ExExFloat ExExFloat::operator =(double value){
    mantissa = value;
    exponent = 0;
    normalize();
    return *this;
}
//==============================================================================
// ARITHMETIC OPERATORS
//==============================================================================

// Addition: align exponents, add mantissas, normalize
ExExFloat ExExFloat::operator +(const ExExFloat& other) const{
    ExExFloat result;

    if(other.exponent >= exponent){
        if(other.mantissa == 0.0){
            result.mantissa = mantissa;
            result.exponent = exponent;
            result.normalize();
        } else{
            // Align exponents and add mantissas
            result.mantissa = other.mantissa + mantissa * invPow2[other.exponent - exponent];
            result.exponent = other.exponent;
            result.normalize();
        }
    } else{
        if(mantissa == 0.0){
            result.mantissa = other.mantissa;
            result.exponent = other.exponent;
            result.normalize();
        } else{
            // Align exponents and add mantissas
            result.mantissa = mantissa + other.mantissa * invPow2[exponent - other.exponent];
            result.exponent = exponent;
            result.normalize();
        }
    }
    return result;
}
// Subtraction: align exponents, subtract mantissas, normalize
ExExFloat ExExFloat::operator -(const ExExFloat& other) const{
    ExExFloat result;

    if(other.exponent >= exponent){
        if(other.mantissa == 0.0){
            result.mantissa = mantissa;
            result.exponent = exponent;
            result.normalize();
        } else{
            // Align exponents and subtract mantissas
            result.mantissa = mantissa * invPow2[other.exponent - exponent] - other.mantissa;
            result.exponent = other.exponent;
            result.normalize();
        }
    } else{
        if(mantissa == 0.0){
            result.mantissa = -other.mantissa;
            result.exponent = other.exponent;
            result.normalize();
        } else{
            // Align exponents and subtract mantissas
            result.mantissa = mantissa - other.mantissa * invPow2[exponent - other.exponent];
            result.exponent = exponent;
            result.normalize();
        }
    }
    return result;
}
// Addition assignment: add and normalize in place
ExExFloat ExExFloat::operator +=(ExExFloat const &obj){
	if(obj.exponent >= exponent){
		if(obj.mantissa == 0.0){
			normalize();
		} else{
			mantissa = obj.mantissa + mantissa*invPow2[obj.exponent - exponent];
			exponent = obj.exponent; normalize();
		}
	} else{
		if(mantissa == 0.0){
			mantissa = obj.mantissa; exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa + obj.mantissa*invPow2[exponent - obj.exponent];
			normalize();
		}
	}
	return *this;
}
// Subtraction assignment: subtract and normalize in place
ExExFloat ExExFloat::operator -=(ExExFloat const &obj){
	if(obj.exponent >= exponent){
		if(obj.mantissa == 0.0){
			normalize();
		} else{
			mantissa = mantissa*invPow2[obj.exponent - exponent] - obj.mantissa;
			exponent = obj.exponent; normalize();
		}
	} else{
		if(mantissa == 0.0){
			mantissa = -obj.mantissa; exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa - obj.mantissa*invPow2[exponent - obj.exponent];
			normalize();
		}
	}
	return *this;
}
// Multiplication: multiply mantissas, add exponents, normalize
ExExFloat ExExFloat::operator *(const ExExFloat& other) const{
    ExExFloat result;
    result.mantissa = mantissa * other.mantissa;
    result.exponent = exponent + other.exponent;
    result.normalize();
    return result;
}
// Division: divide mantissas, subtract exponents, normalize
ExExFloat ExExFloat::operator /(const ExExFloat& other) const{
    ExExFloat result;
    result.mantissa = mantissa / other.mantissa;
    result.exponent = exponent - other.exponent;
    result.normalize();
    return result;
}
// Scalar multiplication: multiply by double, normalize
ExExFloat ExExFloat::operator *(double value) const{
    ExExFloat result;
    result.mantissa = mantissa * value;
    result.exponent = exponent;
    result.normalize();
    return result;
}

// Scalar division: divide by double, normalize
ExExFloat ExExFloat::operator /(double value) const{
    ExExFloat result;
    result.mantissa = mantissa / value;
    result.exponent = exponent;
    result.normalize();
    return result;
}
// Multiplication assignment: multiply and normalize in place
ExExFloat ExExFloat::operator *=(const ExExFloat& other){
    mantissa *= other.mantissa;
    exponent += other.exponent;
    normalize();
    return *this;
}

// Division assignment: divide and normalize in place
ExExFloat ExExFloat::operator /=(const ExExFloat& other){
    mantissa /= other.mantissa;
    exponent -= other.exponent;
    normalize();
    return *this;
}

// Scalar multiplication assignment: multiply by double in place
ExExFloat ExExFloat::operator *=(double value){
    mantissa *= value;
    normalize();
    return *this;
}

// Scalar division assignment: divide by double in place
ExExFloat ExExFloat::operator /=(double value){
    mantissa /= value;
    normalize();
    return *this;
}
// Double * ExExFloat: convert double to ExExFloat, then multiply
ExExFloat operator *(double lhs, const ExExFloat& rhs){
    ExExFloat result;
    result.mantissa = rhs.mantissa * lhs;
    result.exponent = rhs.exponent;
    result.normalize();
    return result;
}

// Double / ExExFloat: convert double to ExExFloat, then divide
ExExFloat operator /(double lhs, const ExExFloat& rhs){
    ExExFloat result;
    result.mantissa = lhs / rhs.mantissa;
    result.exponent = -rhs.exponent;
    result.normalize();
    return result;
}
// Compare with double: convert double to ExExFloat, then compare
int ExExFloat::operator >=(double value) const{
    if(value == 0) return (mantissa >= 0);
    else{
        ExExFloat R;
        R = value;
        if(exponent > R.exponent) return 1;
        else if((exponent == R.exponent) && (mantissa >= R.mantissa)) return 1;
        else return 0;
    }
}

// Compare with ExExFloat: align exponents, compare mantissas
int ExExFloat::operator >=(const ExExFloat& other) const{
    if(other.mantissa == 0) return (mantissa >= 0);
    else{
        if(exponent > other.exponent) return 1;
        else if((exponent == other.exponent) && (mantissa >= other.mantissa)) return 1;
        else return 0;
    }
}
// Convert to standard double precision
double ExExFloat::getDouble() const {
    return ldexp(mantissa, exponent);
}
// Get sign: -1 (negative), 0 (zero), +1 (positive)
int ExExFloat::sign() const {
    return (mantissa > 0.0) - (mantissa < 0.0);
}
// Return absolute value (positive magnitude)
ExExFloat ExExFloat::absoluteValue() const {
    ExExFloat res;
    res.mantissa = fabs(mantissa);
    res.exponent = exponent;
    return res;
}
// Square root: compute sqrt(mantissa) and halve exponent
ExExFloat ExExFloat::squareRoot() const {
    ExExFloat res;
    if (exponent % 2 == 0) {
        res.mantissa = sqrt(mantissa);
        res.exponent = exponent / 2;
    } else {
        res.mantissa = sqrt(2 * mantissa);
        res.exponent = (exponent - 1) / 2;
    }
    res.normalize();
    return res;
}

// Initialize global lookup table: precompute 1/2^i for i=0 to MAX_EXP_RANGE
void divdiff_init() {
    invPow2 = new double[MAX_EXP_RANGE];
    double curr = 1;
    for (int i = 0; i < MAX_EXP_RANGE; i++) {
        invPow2[i] = curr;
        curr /= 2;
    }
}

// Free memory allocated for inverse powers of 2 lookup table
// Clean up global resources
void divdiff_clear_up() noexcept {
    delete[] invPow2;
}

// Calculate arithmetic mean of array values
double divdiff::calculateMean(double* values, int count) {
    double sum = 0;
    int i;
    for (i = 0; i < count; i++) {
        sum += values[i];
    }
    return sum / count;
}

// Find maximum absolute difference between any two values
double divdiff::calculateMaxAbsDiff(double* values, int length) {
    double values_max = values[0], values_min = values[0];
    int i;
    for (i = 1; i < length; i++) {
        values_min = min(values_min, values[i]);
        values_max = max(values_max, values[i]);
    }
    return fabs(values_max - values_min);
}

// Get total memory usage in bytes for all allocated arrays
long long int divdiff::getMemoryUsage() {
    long long int sum = 0;
    sum += sizeof(double) * maxLen;
    sum += 2 * sizeof(ExExFloat) * (maxLen + EXTRA_LEN + 1);
    sum += sizeof(ExExFloat) * maxLen * maxScale;
    sum += 2 * sizeof(double*) + 3 * sizeof(ExExFloat*) + 4 * sizeof(int) + sizeof(double) + sizeof(ExExFloat);
    return sum;
}

// Constructor: initialize with max points and scaling factor
divdiff::divdiff(int maxLen_, int maxScale_) { // constructor
    maxLen = maxLen_;
    maxScale = maxScale_;
    allocateMemory();
    if (invPow2 == NULL) {
        printf("Error: invPow2 has not been initialized\n");
        exit(EXIT_FAILURE);
    }
}

// Copy constructor: deep copy all arrays and state
divdiff::divdiff(const divdiff& other) { // copy constructor
    maxLen = other.maxLen;
    maxScale = other.maxScale;
    allocateMemory();
    currentLength = other.currentLength;
    scaleFactor = other.scaleFactor;
    meanVal = other.meanVal;
    expMu = other.expMu;
    memcpy(zVals, other.zVals, maxLen * sizeof(double));
    memcpy(hCoeffs, other.hCoeffs, (maxLen + EXTRA_LEN + 1) * sizeof(ExExFloat));
    memcpy(divDiffs, other.divDiffs, (maxLen + EXTRA_LEN + 1) * sizeof(ExExFloat));
    memcpy(derivTerms, other.derivTerms, maxLen * maxScale * sizeof(ExExFloat));
}

// Copy assignment: free old memory, then deep copy from other
divdiff& divdiff::operator=(const divdiff& other){ // copy assignment operator
	freeMemory();
	maxLen=other.maxLen; maxScale=other.maxScale; allocateMemory();
	currentLength=other.currentLength; scaleFactor=other.scaleFactor; meanVal=other.meanVal; expMu=other.expMu;
	memcpy(zVals,other.zVals,maxLen*sizeof(double));
	memcpy(hCoeffs,other.hCoeffs,(maxLen+EXTRA_LEN+1)*sizeof(ExExFloat));
	memcpy(divDiffs,other.divDiffs,(maxLen+EXTRA_LEN+1)*sizeof(ExExFloat));
	memcpy(derivTerms,other.derivTerms,maxLen*maxScale*sizeof(ExExFloat));
	return *this;
}

// Destructor: free all allocated memory
divdiff::~divdiff(){ // destructor
	freeMemory();
}

// Allocate memory for all internal arrays and initialize state
void divdiff::allocateMemory(){
	zVals = new double[maxLen]; hCoeffs = new ExExFloat[maxLen+EXTRA_LEN+1];
	divDiffs = new ExExFloat[maxLen+EXTRA_LEN+1]; derivTerms = new ExExFloat[maxLen*maxScale];
	currentLength = 0; scaleFactor = 1;
}

// Free all allocated memory for internal arrays
void divdiff::freeMemory(){
	delete[] zVals; delete[] hCoeffs; delete[] divDiffs; delete[] derivTerms;
}

// Print ExExFloat array for debugging
void divdiff::printExExFloatList(ExExFloat* list, int length, const char* name){
	int i;
	printf("%s={",name);
	for(i=0;i<length;i++){
		list[i].print();
		if(i<length-1) printf(",");
	}
	printf("};\n");
}

// Print double array for debugging
void divdiff::printDoubleList(double* list, int length, const char* name){
	int i;
	printf("%s={",name);
	for(i=0;i<length;i++){
		printf("%.17g",list[i]);
		if(i<length-1) printf(",");
	}
	printf("};\n");
}

// Print integer array for debugging
void divdiff::printIntegerList(int* list, int length, const char* name){
	int i;
	printf("%s={",name);
	for(i=0;i<length;i++){
		printf("%d",list[i]);
		if(i<length-1) printf(",");
	}
	printf("};\n");
}

// Print vector contents for debugging
void divdiff::printVectorList(std::vector<int> list, const char* name){
	int i; int length = list.size();
	printf("%s={",name);
	for(i=0;i<length;i++){
		printf("%d",list[i]);
		if(i<length-1) printf(",");
	}
	printf("};\n");
}

// Backup z array to temporary storage for reallocation
void divdiff::backupZValues(int length) {
    zBackups = new double[length];
    memcpy(zBackups, zVals, length * sizeof(double));
}
// Restore z array from backup after reallocation
void divdiff::restoreZValues(int length) {
    memcpy(zVals, zBackups, length * sizeof(double));
    delete[] zBackups;
}

// Check if scaling factor needs adjustment due to extreme values
int divdiff::scalingChanged() {
    return fabs(zVals[currentLength - 1] - meanVal) / SCAL_THRESHOLD > scaleFactor;
}

// Add new point and update divided differences efficiently
void divdiff::addElement(double newZValue, int forcedScale, double forcedCentral) {
    int j, k, n, N;
    ExExFloat curr;
    n = currentLength;
    N = maxLen + EXTRA_LEN;
    zVals[n] = newZValue;
    currentLength++;

    if (currentLength == 1) {
        scaleFactor = (forcedScale == 0) ? 1 : forcedScale;
        meanVal = (forcedCentral == 0) ? zVals[0] : forcedCentral;
        expMu.initExpMu(meanVal);
        hCoeffs[0] = 1;
        for (k = 1; k <= N; k++) {
            hCoeffs[k] = hCoeffs[k - 1] / scaleFactor;
        }
        if (meanVal != zVals[0]) {
            for (k = N; k > 0; k--) {
                hCoeffs[k - 1] += hCoeffs[k] * (zVals[0] - meanVal) / k;
            }
        }
        curr = expMu * hCoeffs[0];
        for (k = 0; k < scaleFactor - 1; k++) {
            derivTerms[k * maxLen] = curr;
            curr *= hCoeffs[0];
        }
        divDiffs[0].initExpMu(zVals[0]); // alternatively: divDiffs[0] = curr;
    } else if (scalingChanged() || (currentLength >= maxLen)) {
        addAllElements(currentLength, forcedScale);
    } else {
        for (k = N; k > n; k--) {
            hCoeffs[k - 1] += hCoeffs[k] * (zVals[n] - meanVal) / k;
        }
        curr = expMu * hCoeffs[n];
        for (k = n; k >= 1; k--) {
            hCoeffs[k - 1] = (hCoeffs[k - 1] * n + hCoeffs[k] * (zVals[n] - zVals[n - k])) / (n - k + 1);
        }
        for (k = 0; k < scaleFactor - 1; k++) {
            derivTerms[k * maxLen + n] = curr;
            curr = derivTerms[k * maxLen] * hCoeffs[n];
            for (j = 1; j <= n; j++) {
                curr += derivTerms[k * maxLen + j] * hCoeffs[n - j];
            }
        }
        divDiffs[n] = curr;
    }
}

// Remove most recently added element and update calculations
void divdiff::removeElement() {
    int k, n, N;
    if (currentLength >= 1) {
        n = currentLength - 1;
        N = maxLen + EXTRA_LEN;
        for (k = 1; k <= n; k++) {
            hCoeffs[k - 1] = (hCoeffs[k - 1] * (n - k + 1) - hCoeffs[k] * (zVals[n] - zVals[n - k])) / n;
        }
        for (k = n + 1; k <= N; k++) {
            hCoeffs[k - 1] -= hCoeffs[k] * (zVals[n] - meanVal) / k;
        }
        currentLength--;
    }
}

// Remove specific value from bulk (not just last element)
int divdiff::removeValue(double targetValue, int forcedScale, double forcedCentral) { // remove from bulk
    int j, k, n = currentLength - 1, found = 0;
    for (k = n; k >= 0; k--) {
        if (zVals[k] == targetValue) {
            for (j = n; j >= k; j--) {
                removeElement();
            }
            for (j = k; j < n; j++) {
                addElement(zVals[j + 1], forcedScale, forcedCentral);
            }
            found = 1;
            break;
        }
    }
    return found;
}

// Bulk recalculation: recompute all divided differences from input points
void divdiff::addAllElements(int length, int forcedScale) { // input is taken from zVals, output is written to divDiffs, size of zVals should be not smaller than length
    int i, s;
    currentLength = 0;
    if (forcedScale == 0) {
        s = (int)ceil(calculateMaxAbsDiff(zVals, length) / SCAL_THRESHOLD);
    } else {
        s = forcedScale;
    }
    if ((s > maxScale) || (length >= maxLen)) {
        i = maxLen;
        backupZValues(i);
        freeMemory();
        if (s > maxScale) {
            maxScale = max(maxScale * 2, s);
        }
        if (length >= maxLen) {
            maxLen = max(maxLen * 2, length);
        }
        allocateMemory();
        restoreZValues(i);
    }
    addElement(zVals[0], s, calculateMean(zVals, length));
    for (i = 1; i < length; i++) {
        addElement(zVals[i]); // calculates the vector (d[z_0], 1! d[z_0,z_1], ..., n! d[z_0,z_1,...,z_n]).
    }
}