cmake_minimum_required(VERSION 3.16)
project(PMRQMCModern VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/OpSum.cpp
)

# Header files
set(HEADERS
    include/OpSum.h
)

# Example executable
add_executable(pmrqmc examples/main.cpp ${SOURCES})

# Main library (if needed)
add_library(pmrqmc_lib STATIC ${SOURCES} ${HEADERS})
target_include_directories(pmrqmc_lib PUBLIC include)

# Compiler flags
if(MSVC)
    target_compile_options(pmrqmc PRIVATE /W4)
else()
    target_compile_options(pmrqmc PRIVATE -Wall -Wextra -pedantic)
endif()