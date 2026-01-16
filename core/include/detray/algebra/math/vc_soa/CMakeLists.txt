# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_soa_math vc_soa_math
   "include/algebra/math/vc_soa.hpp"
   "include/algebra/math/impl/vc_soa_boolean.hpp"
   "include/algebra/math/impl/vc_soa_math.hpp"
   "include/algebra/math/impl/vc_soa_matrix.hpp"
   "include/algebra/math/impl/vc_soa_vector.hpp")
target_link_libraries( algebra_vc_soa_math
   INTERFACE algebra::common algebra::common_math algebra::common_storage  algebra::vc_soa_storage Vc::Vc )
algebra_test_public_headers( algebra_vc_soa_math
   "algebra/math/vc_soa.hpp" )
