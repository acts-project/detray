/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
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

// Detray test include(s).
#include "detray/test/utils/types.hpp"

namespace detray {

using matadata_t = test::toy_metadata;
using test_algebra = matadata_t::algebra_type;
using scalar = detray::dscalar<test_algebra>;

using detector_host_type =
    detray::detector<matadata_t, detray::host_container_types>;
using detector_device_type =
    detray::detector<matadata_t, detray::device_container_types>;

using navigator_host_type = detray::navigator<detector_host_type>;
using navigator_device_type = detray::navigator<detector_device_type>;
using field_type = detray::bfield::const_field_t<scalar>;
using rk_stepper_type = detray::rk_stepper<field_type::view_t, test_algebra>;
using actor_chain_t =
    detray::actor_chain<detray::parameter_transporter<test_algebra>,
                        detray::pointwise_material_interactor<test_algebra>,
                        detray::parameter_resetter<test_algebra>>;
using propagator_host_type =
    detray::propagator<rk_stepper_type, navigator_host_type, actor_chain_t>;
using propagator_device_type =
    detray::propagator<rk_stepper_type, navigator_device_type, actor_chain_t>;

enum class propagate_option {
    e_unsync = 0,
    e_sync = 1,
};

/// test function for propagator with single state
void propagator_benchmark(
    typename detector_host_type::view_type det_data,
    typename field_type::view_t field_data,
    vecmem::data::vector_view<free_track_parameters<test_algebra>>& tracks_data,
    const propagate_option opt);

}  // namespace detray
