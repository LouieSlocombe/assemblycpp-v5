# tests/CMakeLists.txt

# Gather all test source files
file(GLOB TEST_SOURCES test_*.cpp)
set(TEST_SOURCES ${TEST_SOURCES} main.cpp)

# Create a unit test executable
add_executable(unit_tests ${TEST_SOURCES})

# Link against the core library and Catch2
target_link_libraries(unit_tests
    PRIVATE
        assembly_core
        Catch2::Catch2
)

# Workaround for filesystem linking issues on some GCC setups
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(unit_tests PRIVATE stdc++fs)
endif()

# Include the public headers
target_include_directories(unit_tests PRIVATE ../include)

# Enable test discovery and registration
include(Catch)
catch_discover_tests(unit_tests)
