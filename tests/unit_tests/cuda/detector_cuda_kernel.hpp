/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#if defined(array)
#include "detray/plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#endif

#include <thrust/tuple.h>

#include "detray/core/detector.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

/// test function for detector
void detector_test(detector_data<detector<darray, thrust::tuple, vecmem::vector,
                                          vecmem::jagged_vector> >& det_data);