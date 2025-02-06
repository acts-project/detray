# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

message(STATUS "Building the 'detray::algebra_smatrix' plugin")

# A sanity check.
if(NOT ALGEBRA_PLUGINS_INCLUDE_SMATRIX)
    message(
        WARNING
        "SMatrix not available from Algebra Plugins. "
        "The configuration will likely fail."
    )
endif()

# Set up the library.
detray_add_library( detray_algebra_smatrix algebra_smatrix
   "include/detray/plugins/algebra/smatrix_definitions.hpp"
)
target_link_libraries(
    detray_algebra_smatrix
    INTERFACE algebra::utils algebra::smatrix_smatrix vecmem::core
)
target_compile_definitions(
    detray_algebra_smatrix
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_SMATRIX
)

# Set up tests for the public header(s) of detray::algebra_smatrix.
detray_test_public_headers( detray_algebra_smatrix
   "detray/plugins/algebra/smatrix_definitions.hpp"
)
