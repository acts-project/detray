/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

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

#include "detray/core/surfaces_finder.hpp"

using namespace detray;
using namespace __plugin;

namespace detray {

constexpr int n_grids = 2;

using surfaces_finder_host_t =
    surfaces_finder<n_grids, darray, thrust::tuple, vecmem::vector,
                    vecmem::jagged_vector>;

using surfaces_finder_device_t =
    surfaces_finder<n_grids, darray, thrust::tuple, vecmem::device_vector,
                    vecmem::jagged_device_vector>;

using surfaces_regular_circular_grid_t =
    typename surfaces_finder_host_t::surfaces_regular_circular_grid;

void surfaces_finder_test(
    surfaces_finder_view<surfaces_finder_host_t> finder_data,
    vecmem::data::vector_view<dindex>& outputs_data);

}  // namespace detray
