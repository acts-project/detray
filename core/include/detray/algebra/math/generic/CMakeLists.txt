# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library(algebra_generic_math generic_math
   # impl include
   "include/algebra/math/generic.hpp"
   "include/algebra/math/impl/generic_matrix.hpp"
   "include/algebra/math/impl/generic_transform3.hpp"
   "include/algebra/math/impl/generic_vector.hpp"
   # algorithms include
   "include/algebra/math/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
   "include/algebra/math/algorithms/matrix/determinant/cofactor.hpp"
   "include/algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
   "include/algebra/math/algorithms/matrix/determinant/partial_pivot_lud.hpp"
   "include/algebra/math/algorithms/matrix/inverse/cofactor.hpp"
   "include/algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
   "include/algebra/math/algorithms/matrix/inverse/partial_pivot_lud.hpp"
   "include/algebra/math/algorithms/utils/algorithm_finder.hpp")
target_link_libraries(algebra_generic_math
   INTERFACE algebra::common algebra::utils algebra::common_math)
algebra_test_public_headers( algebra_generic_math
   "algebra/math/generic.hpp" )
