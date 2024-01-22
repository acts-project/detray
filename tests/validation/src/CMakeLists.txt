#Detray library, part of the ACTS project(R& D line)
#
#(c) 2023 CERN for the benefit of the ACTS project
#
#Mozilla Public License Version 2.0

# Silence actsvg internal warnings
include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::meta,INTERFACE_INCLUDE_DIRECTORIES> )
include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::core,INTERFACE_INCLUDE_DIRECTORIES> )

# Build the visualization executable.
detray_add_executable(detector_display
                      "detector_display.cpp"
                      LINK_LIBRARIES GTest::gtest GTest::gtest_main  Boost::program_options
                      detray::utils detray::core_array detray::svgtools)

# Build the validation executable.
detray_add_executable(detector_validation
                      "detector_validation.cpp"
                      LINK_LIBRARIES GTest::gtest GTest::gtest_main 
                      detray::validation)

detray_add_executable(material_validation
                      "material_validation.cpp"
                      LINK_LIBRARIES GTest::gtest GTest::gtest_main 
                      detray::validation)

# Build the CI tests
detray_add_test(jacobian_validation
                "jacobian_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main
                detray::validation)

detray_add_test(telescope_detector
                "telescope_detector_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main
                detray::validation)

detray_add_test(toy_detector
                "toy_detector_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main
                detray::validation)

detray_add_test(wire_chamber
                "wire_chamber_validation.cpp"
                LINK_LIBRARIES GTest::gtest GTest::gtest_main
                detray::validation)
