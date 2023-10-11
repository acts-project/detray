/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tutorial/types.hpp"

// Covfie include(s)
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>

namespace detray::tutorial {

// Detector
using detector_host_t = detector<detray::toy_metadata, host_container_types>;
using detector_device_t =
    detector<detray::toy_metadata, device_container_types>;

// Navigator
using navigator_t = navigator<detector_device_t>;
using intersection_t = navigator_t::intersection_type;

// Magnetic field types
using inhom_cuda_bknd_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::cuda_device_array<
            covfie::vector::vector_d<detray::scalar, 3>>>>>;
using field_host_t = covfie::field<detray::bfield::inhom_bknd_t>;
using field_device_t = covfie::field<inhom_cuda_bknd_t>;

// Stepper (only running device propagation, so need the device field type)
using stepper_t =
    rk_stepper<field_device_t::view_t, detray::tutorial::transform3>;

// Actors
using actor_chain_t =
    actor_chain<tuple, pathlimit_aborter,
                parameter_transporter<detray::tutorial::transform3>,
                pointwise_material_interactor<detray::tutorial::transform3>,
                parameter_resetter<detray::tutorial::transform3>>;

// Propagator
using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

/// Propagation tutorial function
void propagation(
    typename detector_host_t::view_type det_data,
    covfie::field_view<inhom_cuda_bknd_t> field_data,
    const vecmem::data::vector_view<
        free_track_parameters<detray::tutorial::transform3>>
        tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data);

}  // namespace detray::tutorial
