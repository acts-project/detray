# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

if(DETRAY_BUILD_HOST)
    add_subdirectory(cpu)
endif()

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()
