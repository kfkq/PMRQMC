# Define a library for our QMC code.
# For now, it's an INTERFACE library because all our code is header-only.
# This makes its include directories available to any target that links against it.
# We will change this to a regular library once we add .cpp files.
add_library(pmrqmc_core INTERFACE)

# Add the 'core' directory to the public include path for this library.
# This means any executable linking against 'pmrqmc_core' can do #include "core/ExExFloat.hpp"
target_include_directories(pmrqmc_core
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)