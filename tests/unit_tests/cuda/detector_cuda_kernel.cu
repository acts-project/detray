/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detector_cuda_kernel.hpp"
#include "detray/definitions/cuda_defs.hpp"

namespace detray {

void detector_test(detector_data<detector<darray, thrust::tuple, vecmem::vector,
                                          vecmem::jagged_vector> >& det_data) {}

}  // namespace detray