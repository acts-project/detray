# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Temporary setting for the scalar type, until it can be removed.
set(DETRAY_CUSTOM_SCALARTYPE
    "double"
    CACHE STRING
    "Scalar type to use in the Detray code"
)

# Add all subdirectories.
add_subdirectory(array)
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()
if(DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
if(DETRAY_VC_AOS_PLUGIN)
    add_subdirectory(vc_aos)
endif()
if(DETRAY_VC_SOA_PLUGIN)
    add_subdirectory(vc_soa)
endif()
