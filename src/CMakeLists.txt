# Define a library for our QMC code.
add_library(pmrqmc_core INTERFACE)

# Add the 'core' directory to the public include path for this library.
# This means any executable linking against 'pmrqmc_core' can do #include "core/ExExFloat.hpp"
target_include_directories(pmrqmc_core
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Require C++20 standard
target_compile_features(pmrqmc_core INTERFACE cxx_std_20)