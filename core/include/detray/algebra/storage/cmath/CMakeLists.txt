# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_cmath_storage cmath_storage
   "include/algebra/storage/impl/cmath_getter.hpp" )
target_link_libraries( algebra_cmath_storage
   INTERFACE algebra::common algebra::common_math )
algebra_test_public_headers( algebra_cmath_storage
   "algebra/storage/impl/cmath_getter.hpp" )
