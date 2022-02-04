# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include all of the benchmarks that can be built.
add_subdirectory( core )
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
