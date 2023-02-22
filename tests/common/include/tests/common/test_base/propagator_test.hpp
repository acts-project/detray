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
#include "detray/simulation/track_generators.hpp"
#include "detray/tracks/tracks.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/vector.hpp>

// GTest include(s).
#include <gtest/gtest.h>

namespace detray {

// Constant magnetic field
using const_bfield_bknd_t = covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                covfie::vector::vector_d<scalar, 3>>;

// Inhomogeneous magnetic field
using inhom_bfield_bcknd_t = covfie::backend::affine<
        covfie::backend::nearest_neighbour<covfie::backend::strided<
            covfie::vector::ulong3,
            covfie::backend::array<covfie::vector::vector_d<scalar, 3>>>>>;

// Host detector type
template<typename bfield_bknd_t>
using detector_host_t = detector<detector_registry::template toy_detector<bfield_bknd_t>, covfie::field, host_container_types>;

// Device detector type using views
template<typename bfield_bknd_t>
using detector_device_t = detector<detector_registry::template toy_detector<bfield_bknd_t>, covfie::field_view, device_container_types>;

// These types are identical in host and device code for all bfield types
using transform3 = typename detector_host_t<const_bfield_bknd_t>::transform3;
using matrix_operator = standard_matrix_operator<scalar>;
using free_track_parameters_type = free_track_parameters<transform3>;
using free_matrix = typename free_track_parameters_type::covariance_type;

// Navigator
using intersection_t = line_plane_intersection;
template<typename detector_t>
using navigator_t = navigator<detector_t>;

// Stepper
using constraints_t = constrained_step<>;
template<typename detector_t>
using rk_stepper_t =
    rk_stepper<typename detector_t::bfield_type::view_t, transform3, constraints_t>;

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
    actor_chain<thrust::tuple, inspector_host_t, pathlimit_aborter,
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;
using actor_chain_device_t =
    actor_chain<thrust::tuple, inspector_device_t, pathlimit_aborter,
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;

template<typename detector_t, typename actor_chain_t>
using propagator_t = propagator<rk_stepper_t<detector_t>, navigator_t<detector_t>, actor_chain_t>;

}  // namespace detray