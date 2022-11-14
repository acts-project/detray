# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the common C++ flags.
include(detray-compiler-options-cpp)

# Include all of the code-holding sub-directories.
add_subdirectory(common)
add_subdirectory(unit_tests)
add_subdirectory(simulation)

if(DETRAY_BENCHMARKS)
   add_subdirectory(benchmarks)
endif()
