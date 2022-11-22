# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.
detray_add_test( io
   "io_json_payload.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io )