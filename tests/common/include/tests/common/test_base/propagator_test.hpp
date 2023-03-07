/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
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
using const_bfield_bknd_t =
    covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                              covfie::vector::vector_d<scalar, 3>>;

// Inhomogeneous magnetic field
using inhom_bfield_bknd_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::vector_d<scalar, 3>>>>>;

// Host detector type
template <typename bfield_bknd_t>
using detector_host_t =
    detector<detector_registry::template toy_detector<bfield_bknd_t>,
             covfie::field, host_container_types>;

// Device detector type using views
template <typename bfield_bknd_t>
using detector_device_t =
    detector<detector_registry::template toy_detector<bfield_bknd_t>,
             covfie::field_view, device_container_types>;

// These types are identical in host and device code for all bfield types
using transform3 = typename detector_host_t<const_bfield_bknd_t>::transform3;
using vector3_t = typename transform3::vector3;
using matrix_operator = standard_matrix_operator<scalar>;
using track_t = free_track_parameters<transform3>;
using free_matrix = typename track_t::covariance_type;

// Navigator
template <typename detector_t>
using navigator_t = navigator<detector_t>;

// Stepper
using constraints_t = constrained_step<>;
template <typename bfield_view_t>
using rk_stepper_t = rk_stepper<bfield_view_t, transform3, constraints_t>;

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

        state(vecmem::memory_resource& resource)
            : _path_lengths(&resource),
              _positions(&resource),
              _jac_transports(&resource) {}

        DETRAY_HOST_DEVICE
        state(vector_t<scalar> path_lengths, vector_t<vector3_t> positions,
              vector_t<free_matrix> jac_transports)
            : _path_lengths(path_lengths),
              _positions(positions),
              _jac_transports(jac_transports) {}

        vector_t<scalar> _path_lengths;
        vector_t<vector3_t> _positions;
        vector_t<free_matrix> _jac_transports;
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state& inspector_state, const propagator_state_t& prop_state) const {

        const auto& stepping = prop_state._stepping;

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
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;
using actor_chain_device_t =
    actor_chain<tuple, inspector_device_t, pathlimit_aborter,
                parameter_transporter<transform3>,
                pointwise_material_interactor<transform3>,
                parameter_resetter<transform3>>;

/// Precompute the tracks
inline vecmem::vector<track_t> generate_tracks(
    vecmem::memory_resource* mr, const unsigned int ts = theta_steps,
    const unsigned int ps = phi_steps) {

    // Putput collection
    vecmem::vector<track_t> tracks(mr);

    // Set origin position of tracks
    const point3 ori{0.f, 0.f, 0.f};
    const scalar p_mag{10.f * unit<scalar>::GeV};

    // Iterate through uniformly distributed momentum directions
    for (auto track : uniform_track_generator<track_t>(
             ts, ps, ori, p_mag, {0.01f, constant<scalar>::pi},
             {-constant<scalar>::pi, constant<scalar>::pi})) {
        track.set_overstep_tolerance(overstep_tolerance);

        // Put it into vector of trajectories
        tracks.push_back(track);
    }

    return tracks;
}

/// test function for propagator on the host
template <typename bfield_bknd_t>
inline auto run_propagation_host(vecmem::memory_resource* mr,
                                 const detector_host_t<bfield_bknd_t>& det,
                                 const vecmem::vector<track_t>& tracks)
    -> std::tuple<vecmem::jagged_vector<scalar>,
                  vecmem::jagged_vector<vector3_t>,
                  vecmem::jagged_vector<free_matrix>> {

    // Construct propagator from stepper and navigator
    auto stepr = rk_stepper_t<
        typename detector_host_t<bfield_bknd_t>::bfield_type::view_t>{};
    auto nav = navigator_t<detector_host_t<bfield_bknd_t>>{};

    using propagator_host_t =
        propagator<decltype(stepr), decltype(nav), actor_chain_host_t>;
    propagator_host_t p(std::move(stepr), std::move(nav));

    // Create vector for track recording
    vecmem::jagged_vector<scalar> host_path_lengths(mr);
    vecmem::jagged_vector<vector3_t> host_positions(mr);
    vecmem::jagged_vector<free_matrix> host_jac_transports(mr);

    for (const auto& trk : tracks) {

        // Create the propagator state
        inspector_host_t::state insp_state{*mr};
        pathlimit_aborter::state pathlimit_state{path_limit};
        parameter_transporter<transform3>::state transporter_state{};
        pointwise_material_interactor<transform3>::state interactor_state{};
        parameter_resetter<transform3>::state resetter_state{};
        auto actor_states =
            detray::tie(insp_state, pathlimit_state, transporter_state,
                        interactor_state, resetter_state);

        typename propagator_host_t::state state(trk, det.get_bfield(), det);

        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            constrainted_step_size);

        state._stepping.set_tolerance(rk_tolerance);

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

/// test function for propagator on the host
inline void compare_propagation_results(
    const vecmem::jagged_vector<vector3_t>& host_positions,
    const vecmem::jagged_vector<vector3_t>& device_positions,
    const vecmem::jagged_vector<scalar>& host_path_lengths,
    const vecmem::jagged_vector<scalar>& device_path_lengths,
    const vecmem::jagged_vector<free_matrix>& host_jac_transports,
    const vecmem::jagged_vector<free_matrix>& device_jac_transports) {

    // Compare the positions
    for (unsigned int i = 0; i < host_positions.size(); i++) {
        ASSERT_TRUE(host_positions[i].size() > 0);

        for (unsigned int j = 0; j < host_positions[i].size(); j++) {

            auto host_pl = host_path_lengths[i][j];
            auto device_pl = device_path_lengths[i][j];

            // ASSERT_NEAR((host_pl - device_pl) / host_pl, 0, is_close);

            ASSERT_EQ(host_positions[i].size(), device_positions[i].size());
            ASSERT_NEAR(host_pl, device_pl, host_pl * is_close);

            auto& host_pos = host_positions[i][j];
            auto& device_pos = device_positions[i][j];

            auto relative_error =
                static_cast<point3>(1. / host_pl * (host_pos - device_pos));

            ASSERT_NEAR(getter::norm(relative_error), 0, is_close);
        }
    }

    // Compare the jacobian transports
    for (unsigned int i = 0; i < host_jac_transports.size(); i++) {
        for (unsigned int j = 0; j < host_jac_transports[i].size(); j++) {

            auto& host_J = host_jac_transports[i][j];
            auto& device_J = device_jac_transports[i][j];

            auto pl = host_path_lengths[i][j];

            for (std::size_t row = 0; row < e_free_size; row++) {
                for (std::size_t col = 0; col < e_free_size; col++) {

                    auto host_val = matrix_operator().element(host_J, row, col);

                    auto device_val =
                        matrix_operator().element(device_J, row, col);

                    ASSERT_NEAR((host_val - device_val) / pl, 0, is_close);
                }
            }
        }
    }
}

}  // namespace detray
