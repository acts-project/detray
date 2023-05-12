/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/examples/types.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

namespace detray {

// Detector
using detector_host_t = detector<detector_registry::toy_detector, covfie::field,
                                 host_container_types>;
using detector_device_t = detector<detector_registry::toy_detector,
                                   covfie::field_view, device_container_types>;

// Navigator
using intersection_t = intersection2D<typename detector_device_t::surface_type,
                                      detray::example::transform3>;
using navigator_t = navigator<detector_device_t>;

// Stepper
using field_t = detector_host_t::bfield_type;
using stepper_t =
    rk_stepper<field_t::view_t, detray::example::transform3>;

// Actors
using actor_chain_t =
    actor_chain<tuple, pathlimit_aborter,
                parameter_transporter<detray::example::transform3>,
                pointwise_material_interactor<detray::example::transform3>,
                parameter_resetter<detray::example::transform3>>;

// Propagator
using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

/// Propagation example function
void propagation_example(
    typename detector_host_t::detector_view_type det_data,
    const vecmem::data::vector_view<
        free_track_parameters<detray::example::transform3>>
        tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data);

}  // namespace detray
