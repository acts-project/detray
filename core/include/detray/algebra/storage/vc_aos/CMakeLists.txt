# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_aos_storage vc_aos_storage
   "include/algebra/storage/vc_aos.hpp"
   "include/algebra/storage/impl/vc_aos_approximately_equal.hpp"
   "include/algebra/storage/impl/vc_aos_concepts.hpp"
   "include/algebra/storage/impl/vc_aos_getter.hpp" )
target_link_libraries( algebra_vc_aos_storage
   INTERFACE algebra::common algebra::common_storage Vc::Vc algebra::common_math )
algebra_test_public_headers( algebra_vc_aos_storage
   "algebra/storage/vc_aos.hpp" )
