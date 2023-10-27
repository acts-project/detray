# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

include( CMakeFindDependencyMacro )

find_dependency( Boost COMPONENTS program_options )

add_library( detray_validation INTERFACE )
add_library( detray::validation ALIAS detray_validation )

target_include_directories( detray_validation
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )

target_link_libraries( detray_validation
   INTERFACE GTest::gtest GTest::gtest_main Boost::program_options
             vecmem::core covfie::core detray::test detray_tests_common
             detray::io detray::utils detray::core_array )

add_subdirectory( src )
