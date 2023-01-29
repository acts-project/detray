# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable(run_tracer
   "trace_rectangle.cpp"
   LINK_LIBRARIES detray_tests_common detray::vc_array)
