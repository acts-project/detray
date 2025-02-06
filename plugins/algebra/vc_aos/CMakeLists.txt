# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

message(STATUS "Building the 'detray::algebra_vc_aos' plugin")

# A sanity check.
if(NOT ALGEBRA_PLUGINS_INCLUDE_VC)
    message(
        WARNING
        "Vc not available from Algebra Plugins. "
        "The configuration will likely fail."
    )
endif()

# Set up the library.
detray_add_library( detray_algebra_vc_aos algebra_vc_aos
   "include/detray/plugins/algebra/vc_aos_definitions.hpp"
)
target_link_libraries(
    detray_algebra_vc_aos
    INTERFACE algebra::utils algebra::vc_aos vecmem::core
)
target_compile_definitions(
    detray_algebra_vc_aos
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_VC_AOS
)

# Set up tests for the public header(s) of detray::algebra_vc_aos.
detray_test_public_headers( detray_algebra_vc_aos
   "detray/plugins/algebra/vc_aos_definitions.hpp"
)
