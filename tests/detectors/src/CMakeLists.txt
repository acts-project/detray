#Detray library, part of the ACTS project(R& D line)
#
#(c) 2023 CERN for the benefit of the ACTS project
#
#Mozilla Public License Version 2.0

include( CMakeFindDependencyMacro )

find_dependency( Boost COMPONENTS program_options )

detray_add_executable(generate_toy_detector
                      "generate_toy_detector.cpp"
                      LINK_LIBRARIES Boost::program_options
                      detray::io detray::utils detray::core_array)

detray_add_executable(generate_wire_chamber
                      "generate_wire_chamber.cpp"
                      LINK_LIBRARIES Boost::program_options
                      detray::io detray::utils detray::core_array)

detray_add_executable(generate_telescope_detector
                      "generate_telescope_detector.cpp"
                      LINK_LIBRARIES Boost::program_options
                      detray::io detray::utils detray::core_array)
