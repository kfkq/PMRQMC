cmake_minimum_required(VERSION 3.16)

project(PMRQMC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable testing for the entire project.
enable_testing()

# Try to find GTest, but don't fail if not found
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    message(STATUS "GTest not found, downloading and building...")
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    # For older versions of CMake, we need to explicitly set CMAKE_POLICY_DEFAULT_CMP0077
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(src)
add_subdirectory(tests)