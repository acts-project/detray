# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s) that are algebra agnostic.
detray_add_test( cpu
   "tools_actor_chain.cpp"
   LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core )

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
