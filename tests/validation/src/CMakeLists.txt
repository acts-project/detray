#Detray library, part of the ACTS project(R& D line)
#
#(c) 2023 CERN for the benefit of the ACTS project
#
#Mozilla Public License Version 2.0

# Build the validation executable.
detray_add_executable(detector_validation
                      "detector_validation.cpp"
                      LINK_LIBRARIES GTest::gtest GTest::gtest_main 
                      detray_validation)

# Build the CI tests
detray_add_test(telescope_detector
                "telescope_detector_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main detray_validation)

detray_add_test(toy_detector
                "toy_detector_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main detray_validation)

detray_add_test(wire_chamber
                "wire_chamber_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main detray_validation)