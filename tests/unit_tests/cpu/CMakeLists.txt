# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s) that are algebra agnostic.
detray_add_test( cpu
   "core_barcode.cpp"
   "core_container.cpp"
   "tools_actor_chain.cpp"
   LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core )

# Macro setting up the CPU tests for a specific algebra plugin.
macro( detray_add_cpu_test algebra )
   # Build the test executable.
   detray_add_test( cpu_${algebra}
      "tools_bounding_volume.cpp"
      LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core_${algebra}
                     detray::test )
endmacro()

# Build the array tests.
detray_add_cpu_test( array )

# Build the Eigen tests.
if( DETRAY_EIGEN_PLUGIN )
   detray_add_cpu_test( eigen )
endif()

# Build the SMatrix tests.
if( DETRAY_SMATRIX_PLUGIN )
   detray_add_cpu_test( smatrix )
endif()

# Build the Vc tests.
if( DETRAY_VC_PLUGIN )
   detray_add_cpu_test( vc )
endif()

# Build the "old style" tests.
add_subdirectory( array )
if( DETRAY_EIGEN_PLUGIN )
   add_subdirectory( eigen )
endif()
if( DETRAY_SMATRIX_PLUGIN )
   add_subdirectory( smatrix )
endif()
if( DETRAY_VC_PLUGIN )
   add_subdirectory( vc )
endif()
