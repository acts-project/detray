# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message (STATUS "Building detray HIP unit tests")

cmake_minimum_required(VERSION 3.21) # HIP langauge support requires minimum 3.21 
cmake_policy(VERSION 3.21.3...3.27)
find_package(hip REQUIRED)

# Enable HIP as a language.
enable_language(HIP)

# Set the HIP build flags.
include(detray-compiler-options-hip)



# Independent of algebra-plugin
detray_add_unit_test( hip
       "container_hip_kernel.hpp"
       "container_hip.hip"
       "container_hip_kernel.hip"
   LINK_LIBRARIES GTest::gtest_main vecmem::hip detray::core
)

