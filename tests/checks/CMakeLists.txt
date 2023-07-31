# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

add_library( detray_checks INTERFACE )
target_include_directories( detray_checks
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_link_libraries( detray_checks
   INTERFACE GTest::gtest GTest::gtest_main detray::core_array
             detray::test detray_tests_common covfie::core vecmem::core
             detray::io detray::utils )

add_subdirectory( src )
