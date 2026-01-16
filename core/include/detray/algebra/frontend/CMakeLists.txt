# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all enabled libraries.
add_subdirectory( array_cmath )

if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   add_subdirectory( eigen_generic )
   add_subdirectory( eigen_eigen )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_SMATRIX )
   add_subdirectory( smatrix_generic )
   add_subdirectory( smatrix_smatrix )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_VC )
   add_subdirectory( vc_aos )
   add_subdirectory( vc_aos_generic )
   add_subdirectory( vc_soa )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_FASTOR )
   add_subdirectory( fastor_fastor )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_VECMEM )
   add_subdirectory( vecmem_cmath )
endif()
