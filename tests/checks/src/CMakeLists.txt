# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Build the executable.
detray_add_executable( detector_checker
      "checks_main.cpp"
      LINK_LIBRARIES GTest::gtest GTest::gtest_main detray_checks)
