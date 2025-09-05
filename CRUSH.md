# PMRQMC Codebase Guide

## Build Commands
```bash
# Main project build (CMake)
mkdir build && cd build
cmake ..
make

# Quick test build and run (legacy)
cd legacy
g++ -O3 -std=c++11 -o prepare.bin prepare.cpp
./prepare.bin H.txt A.txt A.txt A.txt A.txt A.txt B.txt B.txt B.txt B.txt B.txt A.txt A.txt A.txt A.txt A.txt
g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp
./PMRQMC.bin > single_thread_output.txt
```

## Test Commands
```bash
# Run all tests
ctest

# Run specific tests
./build/tests/exexfloat_test      # Test ExExFloat performance optimizations
./build/tests/physics_test        # Test physics representation (Pauli strings, Hamiltonian)

# Single test run (legacy)
./test_run.sh

# MPI parallel test (5+ cores enables auto-thermalization)
mpicxx -O3 -std=c++11 -o PMRQMC_mpi.bin PMRQMC_mpi.cpp
mpirun -n 5 ./PMRQMC_mpi.bin > therm_test_output.txt
```

## Code Style Guidelines

### General
- C++20 standard with modern features (noexcept, constexpr, namespace)
- Header-only core library structure for numerical components
- Performance-critical: use -O3 optimization and thread-local caching
- Scientific computing focus with numerical precision concerns
- RAII principles and modern C++ best practices

### Naming Conventions
- Classes: PascalCase (`ExExFloat`, `PauliString`, `OperatorTerm`, `Hamiltonian`)
- Functions/Methods: snake_case (`multiply_operators()`, `get_qubits()`, `from_exp()`)
- Private members: trailing underscore (`mantissa_`, `operators_`, `coefficient_`)
- Namespaces: lowercase (`pmrqmc::core`)
- Enum classes: PascalCase (`PauliOp`)

### Formatting
- `#pragma once` for header guards
- Three-slash comments for multi-line documentation
- Space before/after operators and parentheses
- Empty lines between logical sections
- Template syntax: `template<typename T>` not `template <class T>`

### Best Practices
- Use `noexcept` for non-throwing functions
- Prefer `constexpr` for compile-time constants
- Handle special values (NaN, inf, zero) explicitly in math functions
- Use STL algorithms and containers (`std::map`, `std::vector`, `std::complex`)
- RAII principles for resource management
- Thread-local caching for performance-critical operations

### Error Handling
- Return special values (NaN for invalid sqrt) rather than exceptions
- Use assertions for internal consistency checks
- Document preconditions/postconditions in comments
- Numerical tolerance for floating-point comparisons

### Imports
- Include what you use, no forward declarations
- System headers first (`<cmath>`, `<iostream>`, `<limits>`)
- Project headers with relative paths (`"core/ExExFloat.hpp"`)
- Prefer angle brackets `<>` over quotes `""` for system includes

## Core Physics Components

### ExExFloat (Extended Exponent Float)
- Purpose: Avoid overflow/underflow in QMC weight calculations
- Optimizations: Thread-local power-of-2 caching for arithmetic operations
- Features: `sqrt()`, `abs()`, `square()`, extreme value handling
- Performance: 5-8x speedup for common operations via cached lookups

### PauliString
- Purpose: Represent products of Pauli operators on specific qubits
- Features: Order-independent equality, commutation checking, state application
- Operations: Multiplication (with proper Pauli algebra), hashing
- Usage: `PauliString("X0 Z1 Y2")` for X⊗Z⊗Y on qubits 0,1,2

### OperatorTerm  
- Purpose: Single Hamiltonian term (coefficient × PauliString)
- Features: Diagonal/off-diagonal classification, scalar arithmetic
- Usage: `OperatorTerm(2.5, PauliString("X0 Z1"))` for 2.5 X₀Z₁

### Hamiltonian
- Purpose: Complete Hamiltonian with file I/O and analysis
- Features: Parse from files, Hermitian verification, term separation
- File formats: Both modern (`1.0 Z0`) and legacy (`1.0 0 Z`) support
- Analysis: Max/min coefficients, diagonal/off-diagonal term counts

## Testing Strategy

### Unit Tests
- **ExExFloat**: Arithmetic correctness, cache functionality, extreme values
- **PauliString**: Algebra, commutation, state application, hashing
- **OperatorTerm**: Classification, arithmetic, coefficient handling
- **Hamiltonian**: File I/O, property analysis, compatibility

### Performance Tests
- Cache effectiveness for small exponent differences
- Large coefficient magnitude handling
- File parsing speed for legacy and modern formats

### Integration Tests
- End-to-end Hamiltonian construction and manipulation
- Mixed arithmetic operations across all types
- File format conversion and compatibility