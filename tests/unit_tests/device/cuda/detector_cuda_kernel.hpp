/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Projetc include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

namespace detray {

// some useful type declarations
using metadata_t = test::toy_metadata;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using det_volume_t = typename detector_host_t::volume_type;
using det_surface_t = typename detector_host_t::surface_type;
using transform_t = typename detector_host_t::transform3_type;
using mask_defs = typename detector_host_t::masks;

constexpr auto rectangle_id = mask_defs::id::e_rectangle2;
constexpr auto disc_id = mask_defs::id::e_portal_ring2;
constexpr auto cylinder_id = mask_defs::id::e_portal_cylinder2;

using rectangle_t = typename mask_defs::template get_type<rectangle_id>::type;
using disc_t = typename mask_defs::template get_type<disc_id>::type;
using cylinder_t = typename mask_defs::template get_type<cylinder_id>::type;

/// declaration of a test function for detector
void detector_test(typename detector_host_t::view_type det_data,
                   vecmem::data::vector_view<det_volume_t> volumes_data,
                   vecmem::data::vector_view<det_surface_t> surfaces_data,
                   vecmem::data::vector_view<transform_t> transforms_data,
                   vecmem::data::vector_view<rectangle_t> rectangles_data,
                   vecmem::data::vector_view<disc_t> discs_data,
                   vecmem::data::vector_view<cylinder_t> cylinders_data);

/// declaration of an alignment test function for detector
void detector_alignment_test(
    typename detector_host_t::view_type det_data_static,
    typename detector_host_t::view_type det_data_aligned,
    vecmem::data::vector_view<transform_t> surfacexf_data_static,
    vecmem::data::vector_view<transform_t> surfacexf_data_aligned);

}  // namespace detray
