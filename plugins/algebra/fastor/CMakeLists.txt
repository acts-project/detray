# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the 'detray::algebra_fastor' plugin")

# A sanity check.
if(NOT ALGEBRA_PLUGINS_INCLUDE_FASTOR)
    message(
        WARNING
        "Fastor not available from Algebra Plugins. "
        "The configuration will likely fail."
    )
endif()

# Set up the library.
detray_add_library( detray_algebra_fastor algebra_fastor
   "include/detray/plugins/algebra/fastor_definitions.hpp"
)
target_link_libraries(
    detray_algebra_fastor
    INTERFACE algebra::utils algebra::fastor_fastor vecmem::core
)
target_compile_definitions(
    detray_algebra_fastor
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_FASTOR
)

# Set up tests for the public header(s) of detray::algebra_fastor.
detray_test_public_headers( detray_algebra_fastor
   "detray/plugins/algebra/fastor_definitions.hpp"
)
