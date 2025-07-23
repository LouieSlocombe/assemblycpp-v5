cmake_minimum_required(VERSION 3.10)
project(assembly)

# Enable LTO for all targets if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported)
  message(STATUS "LTO is supported: enabling interprocedural optimization")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(WARNING "LTO not supported: ${ipo_error}")
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Defaulting to Release.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add include directories
include_directories(include)
include_directories(third_party)

# === Build core library ===
set(CORE_SOURCES
    src/assemblyCalculator.cpp
    src/globalPrimitives.cpp
    src/assemblyState.cpp
    src/dagEnumeration.cpp
    src/duplicateMatching.cpp
    src/fragmentation.cpp
    src/graphHashes.cpp
    src/graphio.cpp
    src/help.cpp
    src/improvedBnB.cpp
    src/ioflag.cpp
    src/molfileParser.cpp
    src/molGraph.cpp
    src/pathwayGenerator.cpp
    src/signalHandler.cpp
    src/treeCanon.cpp
    src/vf2.cpp
)

add_library(assembly_core STATIC ${CORE_SOURCES})
target_include_directories(assembly_core PUBLIC include)

# === Build main program ===
add_executable(assembly src/main.cpp)
target_link_libraries(assembly PRIVATE assembly_core)

# === Set up Catch2 for unit tests ===
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.8.1
)
FetchContent_MakeAvailable(Catch2)

enable_testing()
add_subdirectory(tests)
