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

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/masks/masks.hpp"
#include "vecmem/utils/cuda/copy.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

#include "vecmem/containers/device_vector.hpp"

// Thrust include(s)
#include <thrust/tuple.h>

using namespace detray;

using point2 = __plugin::point2<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<scalar>;
const int n_points = 1000;

namespace detray {

using annulus = mask<annulus2D<>>;
using cylinder = mask<cylinder2D<>>;
using rectangle = mask<rectangle2D<>>;
using ring = mask<ring2D<>>;
using single = mask<single3D<>>;
using trapezoid = mask<trapezoid2D<>>;

/** Enumerate different mask types for convenience
 **/
enum mask_ids : unsigned int {
    e_rectangle2 = 0u,
    e_trapezoid2 = 1u,
    e_ring2 = 2u,
    e_cylinder2 = 3u,
    e_single3 = 4u,
    e_annulus2 = 5u,
};

using host_store_type =
    regular_multi_store<mask_ids, empty_context, dtuple, vecmem::vector,
                        rectangle, trapezoid, ring, cylinder, single, annulus>;

using device_store_type =
    regular_multi_store<mask_ids, empty_context, thrust::tuple,
                        vecmem::device_vector, rectangle, trapezoid, ring,
                        cylinder, single, annulus>;

/// test function for mask store
void mask_test(
    typename host_store_type::view_type store_data,
    vecmem::data::vector_view<point2> input_point2_data,
    vecmem::data::jagged_vector_view<intersection::status> output_data);

}  // namespace detray