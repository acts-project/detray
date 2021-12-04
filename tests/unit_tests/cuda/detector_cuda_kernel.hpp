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

namespace detray {

// some useful type declarations
using detector_t =
    detector<darray, thrust::tuple, vecmem::vector, vecmem::jagged_vector>;
using detector_device_t = detector<darray, thrust::tuple, vecmem::device_vector,
                                   vecmem::jagged_device_vector>;
using volume_t = typename detector_t::volume_type;
using surface_t = typename detector_t::surface_type;
using transform_store_t = typename detector_t::transform_store;
using rectangle_t = typename detector_t::rectangle;
using disc_t = typename detector_t::disc;
using cylinder_t = typename detector_t::cylinder;

/// declaration of a test function for detector
void detector_test(
    detector_data<detector_t>& det_data,
    vecmem::data::vector_view<volume_t>& volumes_data,
    vecmem::data::vector_view<surface_t>& surfaces_data,
    static_transform_store_data<transform_store_t>& transforms_data,
    vecmem::data::vector_view<rectangle_t>& rectangles_data,
    vecmem::data::vector_view<disc_t>& discs_data,
    vecmem::data::vector_view<cylinder_t>& cylinders_data);

}  // namespace detray