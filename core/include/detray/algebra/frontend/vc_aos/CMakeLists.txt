# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_aos vc_aos
   "include/algebra/vc_aos.hpp" )
target_link_libraries( algebra_vc_aos
   INTERFACE algebra::common algebra::vc_aos_storage algebra::vc_aos_math algebra::generic_math )
algebra_test_public_headers( algebra_vc_aos
   "algebra/vc_aos.hpp" )
