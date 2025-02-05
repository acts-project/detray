# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set up all of the "host" tests.
if(DETRAY_BUILD_HOST)
    message(STATUS "Building detray host integration tests")

    add_subdirectory(cpu)
    add_subdirectory(io)
endif()

# Set up all of the "device" tests.
if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL)
    add_subdirectory(device)
endif()
