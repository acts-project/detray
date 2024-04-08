/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <tuple>

namespace detray {

// These types are identical in host and device code for all bfield types
using algebra_t = ALGEBRA_PLUGIN<detray::scalar>;
using scalar_t = dscalar<algebra_t>;
using vector3_t = dvector3D<algebra_t>;
using point3_t = dpoint3D<algebra_t>;
using matrix_operator = dmatrix_operator<algebra_t>;
using track_t = free_track_parameters<algebra_t>;
using free_matrix_t = free_matrix<algebra_t>;

// Navigator
template <typename detector_t>
using navigator_t = navigator<detector_t>;
template <typename detector_t>
using intersection_t = typename navigator_t<detector_t>::intersection_type;

// Stepper
using constraints_t = constrained_step<>;
template <typename bfield_view_t>
using rk_stepper_t = rk_stepper<bfield_view_t, algebra_t, constraints_t>;

// Geomery navigation configurations
constexpr unsigned int theta_steps{10u};
constexpr unsigned int phi_steps{10u};

constexpr scalar_t rk_tolerance{1e-4f};
constexpr scalar_t overstep_tolerance{-3.f * unit<scalar_t>::um};
constexpr scalar_t constrainted_step_size{2.f * unit<scalar_t>::mm};
constexpr scalar_t is_close{1e-4f};
constexpr scalar path_limit{5.f * unit<scalar>::m};

template <template <typename...> class vector_t>
struct track_inspector : actor {

    struct state {

        state(vecmem::memory_resource &resource)
            : _path_lengths(&resource),
              _positions(&resource),
              _jac_transports(&resource) {}

        DETRAY_HOST_DEVICE
        state(vector_t<scalar> path_lengths, vector_t<vector3_t> positions,
              vector_t<free_matrix_t> jac_transports)
            : _path_lengths(path_lengths),
              _positions(positions),
              _jac_transports(jac_transports) {}

        vector_t<scalar> _path_lengths;
        vector_t<vector3_t> _positions;
        vector_t<free_matrix_t> _jac_transports;
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

// Assemble actor chain type
using inspector_host_t = track_inspector<vecmem::vector>;
using inspector_device_t = track_inspector<vecmem::device_vector>;
using actor_chain_host_t =
    actor_chain<tuple, inspector_host_t, pathlimit_aborter,
                parameter_transporter<algebra_t>,
                pointwise_material_interactor<algebra_t>,
                parameter_resetter<algebra_t>>;
using actor_chain_device_t =
    actor_chain<tuple, inspector_device_t, pathlimit_aborter,
                parameter_transporter<algebra_t>,
                pointwise_material_interactor<algebra_t>,
                parameter_resetter<algebra_t>>;

/// Precompute the tracks
template <typename track_generator_t = uniform_track_generator<track_t>>
inline auto generate_tracks(
    vecmem::memory_resource *mr,
    const typename track_generator_t::configuration &cfg = {}) {

    // Track collection
    vecmem::vector<typename track_generator_t::track_type> tracks(mr);

    // Iterate through uniformly distributed momentum directions
    for (auto track : track_generator_t{cfg}) {
        // Put it into vector of trajectories
        tracks.push_back(track);
    }

    return tracks;
}

/// Test function for propagator on the host
template <typename bfield_bknd_t, typename host_detector_t>
inline auto run_propagation_host(vecmem::memory_resource *mr,
                                 const host_detector_t &det,
                                 covfie::field<bfield_bknd_t> &field,
                                 const vecmem::vector<track_t> &tracks)
    -> std::tuple<vecmem::jagged_vector<scalar_t>,
                  vecmem::jagged_vector<vector3_t>,
                  vecmem::jagged_vector<free_matrix_t>> {

    // Construct propagator from stepper and navigator
    auto stepr = rk_stepper_t<typename covfie::field<bfield_bknd_t>::view_t>{};
    auto nav = navigator_t<host_detector_t>{};

    using propagator_host_t =
        propagator<decltype(stepr), decltype(nav), actor_chain_host_t>;
    propagation::config<scalar_t> cfg{};
    cfg.navigation.search_window = {3u, 3u};
    cfg.stepping.rk_error_tol = rk_tolerance;
    propagator_host_t p{cfg};

    // Create vector for track recording
    vecmem::jagged_vector<scalar_t> host_path_lengths(mr);
    vecmem::jagged_vector<vector3_t> host_positions(mr);
    vecmem::jagged_vector<free_matrix_t> host_jac_transports(mr);

    for (const auto &trk : tracks) {

        // Create the propagator state
        inspector_host_t::state insp_state{*mr};
        pathlimit_aborter::state pathlimit_state{path_limit};
        parameter_transporter<algebra_t>::state transporter_state{};
        pointwise_material_interactor<algebra_t>::state interactor_state{};
        parameter_resetter<algebra_t>::state resetter_state{};
        auto actor_states =
            detray::tie(insp_state, pathlimit_state, transporter_state,
                        interactor_state, resetter_state);

        typename propagator_host_t::state state(trk, field, det);

        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            constrainted_step_size);

        // Run propagation
        p.propagate(state, actor_states);

        // Record the step information
        host_path_lengths.push_back(insp_state._path_lengths);
        host_positions.push_back(insp_state._positions);
        host_jac_transports.push_back(insp_state._jac_transports);
    }

    return std::make_tuple(std::move(host_path_lengths),
                           std::move(host_positions),
                           std::move(host_jac_transports));
}

/// Compare the results between host and device propagation
inline void compare_propagation_results(
    const vecmem::jagged_vector<vector3_t> &host_positions,
    const vecmem::jagged_vector<vector3_t> &device_positions,
    const vecmem::jagged_vector<scalar_t> &host_path_lengths,
    const vecmem::jagged_vector<scalar_t> &device_path_lengths,
    const vecmem::jagged_vector<free_matrix_t> &host_jac_transports,
    const vecmem::jagged_vector<free_matrix_t> &device_jac_transports) {

    // Compare the positions and pathlengths
    for (unsigned int i = 0u; i < host_positions.size(); i++) {
        ASSERT_TRUE(host_positions[i].size() > 0);

        for (unsigned int j = 0u; j < host_positions[i].size(); j++) {

            scalar_t host_pl = host_path_lengths[i][j];
            scalar_t device_pl = device_path_lengths[i][j];

            ASSERT_EQ(host_positions[i].size(), device_positions[i].size());
            ASSERT_NEAR(host_pl, device_pl, host_pl * is_close);

            const vector3_t &host_pos = host_positions[i][j];
            const vector3_t &device_pos = device_positions[i][j];

            auto relative_error =
                static_cast<point3_t>(1. / host_pl * (host_pos - device_pos));

            ASSERT_NEAR(getter::norm(relative_error), 0.f, is_close);
        }
    }

    // Compare the Jacobians
    for (unsigned int i = 0u; i < host_jac_transports.size(); i++) {
        for (unsigned int j = 0u; j < host_jac_transports[i].size(); j++) {

            const free_matrix_t &host_J = host_jac_transports[i][j];
            const free_matrix_t &device_J = device_jac_transports[i][j];

            scalar_t pl = host_path_lengths[i][j];

            for (std::size_t row = 0u; row < e_free_size; row++) {
                for (std::size_t col = 0u; col < e_free_size; col++) {

                    scalar_t host_val =
                        matrix_operator().element(host_J, row, col);

                    scalar_t device_val =
                        matrix_operator().element(device_J, row, col);

                    ASSERT_NEAR((host_val - device_val) / pl, 0.f, is_close);
                }
            }
        }
    }
}

}  // namespace detray
