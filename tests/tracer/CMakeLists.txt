# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable(run_tracer
   "single_shapes.cpp"
   LINK_LIBRARIES detray_algebra_vc detray::test detray::io detray::utils)
