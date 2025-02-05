# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set up a device test library
add_library(detray_test_device INTERFACE "propagator_test.hpp")
add_library(detray::test_device ALIAS detray_test_device)

target_link_libraries(
    detray_test_device
    INTERFACE detray::core_array detray::test_common detray::test_utils
)

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()
