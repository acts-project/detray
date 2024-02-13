/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

using namespace detray;

using transform3 = __plugin::transform3<scalar>;

using detector_host_type = detector<toy_metadata, host_container_types>;
using detector_device_type = detector<toy_metadata, device_container_types>;

using intersection_t =
    intersection2D<typename detector_device_type::surface_type, transform3>;

using navigator_host_type = navigator<detector_host_type>;
using navigator_device_type = navigator<detector_device_type>;
using field_type = bfield::const_field_t;
using rk_stepper_type = rk_stepper<field_type::view_t, transform3>;
using actor_chain_t = actor_chain<tuple, parameter_transporter<transform3>,
                                  pointwise_material_interactor<transform3>,
                                  parameter_resetter<transform3>>;
using propagator_host_type =
    propagator<rk_stepper_type, navigator_host_type, actor_chain_t>;
using propagator_device_type =
    propagator<rk_stepper_type, navigator_device_type, actor_chain_t>;

enum class propagate_option {
    e_unsync = 0,
    e_sync = 1,
};

namespace detray {

/// test function for propagator with single state
void propagator_benchmark(
    typename detector_host_type::view_type det_data,
    typename field_type::view_t field_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    const propagate_option opt);

}  // namespace detray
