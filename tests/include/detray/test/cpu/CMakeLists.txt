# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set the CUDA build flags.
include(detray-compiler-options-cpp)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    detray_test_cpu
    INTERFACE
    "detector_consistency.hpp"
    "detector_scan.hpp"
    "material_scan.hpp"
    "navigation_validation.hpp"
    "navigation_validation.hpp"
    "toy_detector_test.hpp"
)

add_library(detray::test_cpu ALIAS detray_test_cpu)

target_link_libraries(
    detray_test_cpu
    INTERFACE
        vecmem::core
        detray::core_array
        detray::test_common
        detray::validation_utils
)
