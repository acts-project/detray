/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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

#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

#include "detray/definitions/units.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/navigation_policies.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

using transform3 = __plugin::transform3<scalar>;
using intersection_t = line_plane_intersection;

using detector_host_type =
    detector<detector_registry::toy_detector, covfie::field, darray,
             thrust::tuple, vecmem::vector, vecmem::jagged_vector>;
using detector_device_type =
    detector<detector_registry::toy_detector, covfie::field_view, darray,
             thrust::tuple, vecmem::device_vector,
             vecmem::jagged_device_vector>;

using navigator_host_type = navigator<detector_host_type>;
using navigator_device_type = navigator<detector_device_type>;
using field_type = detector_host_type::bfield_type;
template <typename policy_t>
using rk_stepper_type =
    rk_stepper<field_type::view_t, transform3, unconstrained_step, policy_t>;
using actor_chain_t =
    actor_chain<thrust::tuple, parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;
template <typename stepper_policy_t>
using propagator_host_type = propagator<rk_stepper_type<stepper_policy_t>,
                                        navigator_host_type, actor_chain_t>;
template <typename stepper_policy_t>
using propagator_device_type = propagator<rk_stepper_type<stepper_policy_t>,
                                          navigator_device_type, actor_chain_t>;

namespace detray {

/// test function for propagator with single state
template <typename stepper_policy_t>
void propagator_benchmark(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data);

}  // namespace detray
