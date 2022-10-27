# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable(run_simulation
   "run_simulation.cpp"
   LINK_LIBRARIES detray_tests_common detray::array)
