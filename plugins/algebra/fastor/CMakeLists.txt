# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the Fastor algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_fastor algebra_fastor
   "include/algebra/fastor.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
)

target_link_libraries(
    detray_algebra_fastor
    INTERFACE Fastor::Fastor detray::core
)

target_compile_definitions(
    detray_algebra_fastor
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_FASTOR
)

detray_test_public_headers( detray_algebra_fastor
   "algebra/fastor.hpp"
)
