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
#include "detray/detectors/detector_metadata.hpp"
#include "detray/utils/ranges.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

// some useful type declarations
using detector_host_t = detector<detector_registry::toy_detector, covfie::field,
                                 host_container_types>;
using detector_device_t = detector<detector_registry::toy_detector,
                                   covfie::field_view, device_container_types>;
using volume_t = typename detector_host_t::volume_type;
using surface_t = typename detector_host_t::surface_type;
using transform_t = typename detector_host_t::transform3;
using transform_store_t = typename detector_host_t::transform_container;
using mask_defs = typename detector_host_t::masks;

constexpr auto rectangle_id = mask_defs::id::e_rectangle2;
constexpr auto disc_id = mask_defs::id::e_portal_ring2;
constexpr auto cylinder_id = mask_defs::id::e_portal_cylinder2;

using rectangle_t = typename mask_defs::template get_type<rectangle_id>::type;
using disc_t = typename mask_defs::template get_type<disc_id>::type;
using cylinder_t = typename mask_defs::template get_type<cylinder_id>::type;

/// declaration of a test function for detector
void detector_test(detector_view<detector_host_t> det_data,
                   vecmem::data::vector_view<volume_t> volumes_data,
                   vecmem::data::vector_view<surface_t> surfaces_data,
                   vecmem::data::vector_view<transform_t> transforms_data,
                   vecmem::data::vector_view<rectangle_t> rectangles_data,
                   vecmem::data::vector_view<disc_t> discs_data,
                   vecmem::data::vector_view<cylinder_t> cylinders_data);

// declaration of a test function for volume enumeration
void enumerate_test(detector_view<detector_host_t> det_data,
                    vecmem::data::jagged_vector_view<surface_t> surfaces_data);

}  // namespace detray
