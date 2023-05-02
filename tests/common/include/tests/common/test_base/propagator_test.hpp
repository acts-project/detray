/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
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

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// GTest include(s).
#include <gtest/gtest.h>

namespace detray {

// Type declarations
using transform3 = __plugin::transform3<scalar>;
using detector_host_type = detector<detector_registry::toy_detector,
                                    covfie::field, host_container_types>;
using detector_device_type =
    detector<detector_registry::toy_detector, covfie::field_view,
             device_container_types>;

using intersection_t =
    intersection2D<typename detector_device_type::surface_type, transform3>;

using navigator_host_type = navigator<detector_host_type>;
using navigator_device_type = navigator<detector_device_type>;

using constraints_t = constrained_step<>;
using field_type = detector_host_type::bfield_type;
using rk_stepper_type =
    rk_stepper<field_type::view_t, transform3, constraints_t>;

using matrix_operator = standard_matrix_operator<scalar>;
using free_track_parameters_type = free_track_parameters<transform3>;
using free_matrix = typename free_track_parameters_type::covariance_type;

// Detector configuration
constexpr std::size_t n_brl_layers{4u};
constexpr std::size_t n_edc_layers{3u};

// Geomery navigation configurations
constexpr unsigned int theta_steps{10u};
constexpr unsigned int phi_steps{10u};

constexpr scalar rk_tolerance{1e-4f};
constexpr scalar overstep_tolerance{-3.f * unit<scalar>::um};
constexpr scalar constrainted_step_size{2.f * unit<scalar>::mm};
constexpr scalar is_close{1e-4f};
constexpr scalar path_limit{2.f * unit<scalar>::m};

template <template <typename...> class vector_t>
struct track_inspector : actor {

    struct state {

        state(vecmem::memory_resource &resource)
            : _path_lengths(&resource),
              _positions(&resource),
              _jac_transports(&resource) {}

        DETRAY_HOST_DEVICE
        state(vector_t<scalar> path_lengths, vector_t<vector3> positions,
              vector_t<free_matrix> jac_transports)
            : _path_lengths(path_lengths),
              _positions(positions),
              _jac_transports(jac_transports) {}

        vector_t<scalar> _path_lengths;
        vector_t<vector3> _positions;
        vector_t<free_matrix> _jac_transports;
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state &inspector_state, const propagator_state_t &prop_state) const {

        const auto &stepping = prop_state._stepping;

        // Nothing happened yet: First call of actor chain
        if (stepping.path_length() < is_close) {
            return;
        }

        // Record only on the object
        inspector_state._path_lengths.push_back(stepping.path_length());
        inspector_state._positions.push_back(stepping().pos());
        inspector_state._jac_transports.push_back(stepping._jac_transport);
    }
};

// Assemble propagator type
using inspector_host_t = track_inspector<vecmem::vector>;
using inspector_device_t = track_inspector<vecmem::device_vector>;
using actor_chain_host_t =
    actor_chain<tuple, inspector_host_t, pathlimit_aborter,
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;
using actor_chain_device_t =
    actor_chain<tuple, inspector_device_t, pathlimit_aborter,
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;
using propagator_host_type =
    propagator<rk_stepper_type, navigator_host_type, actor_chain_host_t>;
using propagator_device_type =
    propagator<rk_stepper_type, navigator_device_type, actor_chain_device_t>;

}  // namespace detray
