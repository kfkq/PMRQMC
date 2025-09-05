# Enable testing for ExExFloat and other core components

# Test executable for ExExFloat
add_executable(exexfloat_test
    exexfloat_test.cpp
)

# Link with core library only (no GTest for now)
target_link_libraries(exexfloat_test PRIVATE pmrqmc_core)

# Add test to CTest
add_test(NAME ExExFloatTest COMMAND exexfloat_test)

# Enable C++20 for test
target_compile_features(exexfloat_test PRIVATE cxx_std_20)

# Test executable for Physics Representation
add_executable(physics_test
    physics_test.cpp
)

# Link with core library
target_link_libraries(physics_test PRIVATE pmrqmc_core)

# Add test to CTest
add_test(NAME PhysicsTest COMMAND physics_test)

# Enable C++20 for test
target_compile_features(physics_test PRIVATE cxx_std_20)

# Filesystem library support (C++17+)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.9)
    target_link_libraries(physics_test PRIVATE stdc++fs)
endif()

# Debug executable for PauliString
add_executable(debug_pauli
    debug_pauli.cpp
)

# Link with core library
target_link_libraries(debug_pauli PRIVATE pmrqmc_core)

# Enable C++20 for test
target_compile_features(debug_pauli PRIVATE cxx_std_20)