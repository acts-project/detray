# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set up all of the "host" tutorials.
if(DETRAY_BUILD_HOST)
    add_subdirectory(cpu)
endif()

# Set up all of the "device" tutorials.
if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL)
    add_subdirectory(device)
endif()
