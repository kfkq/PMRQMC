# PMRQMC Project Information

## Build Commands
- **C++ single-threaded build**: `g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp`
- **C++ MPI build**: `mpicxx -O3 -std=c++11 -o PMRQMC_mpi.bin PMRQMC_mpi.cpp`
- **Prepare data**: `g++ -O3 -std=c++11 -o prepare.bin prepare.cpp`
- **Quick test run**: `./test_run.sh`
- **Python experiment drivers**:
  - `python experiments/tfim_driver.py`
  - `python experiments/xxz_driver.py`

## Test Commands
- **Python unit tests**: `python -m unittest utils/_pauli_manipulations_test.py`
- **Single test method**: `python -m unittest utils._pauli_manipulations_test.TestClass.test_method`
- **Run basic simulation**: `./test_run.sh` (generates single_thread_output.txt)

No dedicated test runner or CI setup. Tests are primarily Python unittest-based.

## Code Style Guidelines

### C++ Conventions
- **Standard**: C++11 (`-std=c++11`)
- **Optimization**: Use `-O3` for production builds
- **Naming**: Use descriptive variable names, consistent with scientific computing style
- **Error Handling**: Use `std::cout` for errors, `exit(1)` on failure
- **Includes**: Local includes use quotes `""`, system headers use angle brackets `<>`
- **Comments**: C-style comments with detailed paper references and licensing info
- **Formatting**: Mixed tabs/spaces usage observed, favor consistency within files

### Python Conventions
- **Imports**: 
  - Standard library first, then third-party (numpy), then local modules
  - Use `sys.path.append("../utils")` for relative imports when needed
  - Avoid wildcard imports except for specific cases
- **Naming**: snake_case for variables/functions, CamelCase for classes
- **Docstrings**: Triple-quoted strings for module/class/function documentation
- **Error Handling**: Standard try/except patterns, raise appropriate exceptions
- **Code Organization**: Keep utility functions in dedicated modules
- **Style**: Follow PEP 8, use type hints when beneficial for clarity

### General Guidelines
- **Documentation**: Include paper references in file headers
- **Modularity**: Separate concerns (C++ for QMC, Python for drivers/utilities)
- **Data Handling**: Use text files for configuration (H.txt, A.txt, B.txt, parameters.hpp)
- **Version Control**: Commit working configurations rather than generated parameters.hpp
- **Reproducibility**: Seed random number generators appropriately for reproducible results</content>
<parameter name="file_path">CRUSH.md