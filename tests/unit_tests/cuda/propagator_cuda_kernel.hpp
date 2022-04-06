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

#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

using intersection_t = line_plane_intersection;

using detector_host_type =
    detector<detector_registry::toy_detector, darray, thrust::tuple,
             vecmem::vector, vecmem::jagged_vector>;
using detector_device_type =
    detector<detector_registry::toy_detector, darray, thrust::tuple,
             vecmem::device_vector, vecmem::jagged_device_vector>;

using navigator_host_type = navigator<detector_host_type>;
using navigator_device_type = navigator<detector_device_type>;

using field_type = constant_magnetic_field<>;
using rk_stepper_type = rk_stepper<field_type, free_track_parameters>;

// detector configuration
constexpr std::size_t n_brl_layers = 4;
constexpr std::size_t n_edc_layers = 3;

// geomery navigation configurations
constexpr unsigned int theta_steps = 100;
constexpr unsigned int phi_steps = 100;

constexpr scalar pos_diff_tolerance = 1e-3;

namespace detray {

template <template <typename...> class vector_t>
struct track_inspector : actor {

    struct track_inspector_state {

        track_inspector_state(vecmem::memory_resource &resource)
            : _intersections(&resource) {}

        DETRAY_HOST_DEVICE
        track_inspector_state(vector_t<intersection_t> intersection_record)
            : _intersections(intersection_record) {}

        vector_t<intersection_t> _intersections;
    };

    using state_type = track_inspector_state;

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state_type &inspector_state,
        const propagator_state_t &prop_state) const {

        const auto &navigation = prop_state._navigation;

        // Record when status == e_on_target
        if (navigation.status() == navigation::status::e_on_target) {
            inspector_state._intersections.push_back(*navigation.current());
        }
    }
};

// Assemble propagator type
using inspector_host_t = track_inspector<vecmem::vector>;
using inspector_device_t = track_inspector<vecmem::device_vector>;
using actor_chain_host_t = actor_chain<thrust::tuple, inspector_host_t>;
using actor_chain_device_t = actor_chain<thrust::tuple, inspector_device_t>;
using propagator_host_type =
    propagator<rk_stepper_type, navigator_host_type, actor_chain_host_t>;
using propagator_device_type =
    propagator<rk_stepper_type, navigator_device_type, actor_chain_device_t>;

/// test function for propagator with single state
void propagator_test(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters> &tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> &candidates_data,
    vecmem::data::jagged_vector_view<intersection_t> &intersections_data);

}  // namespace detray
