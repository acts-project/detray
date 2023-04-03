/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s).
#include <random>

namespace detray {

template <typename interactor_t>
struct random_scatterer : actor {
    using interactor_type = interactor_t;
    using interactor_state_type = typename interactor_type::state;
    using matrix_operator = typename interactor_type::matrix_operator;
    using scalar_type = typename interactor_type::scalar_type;
    using vector3 = typename interactor_type::vector3;
    using transform3_type = typename interactor_type::transform3_type;

    struct state {
        std::random_device rd{};
        std::mt19937_64 generator{rd()};

        /// Constructor with seed
        ///
        /// @param sd the seed number
        state(const uint_fast64_t sd = 0u) { generator.seed(sd); }

        void set_seed(const uint_fast64_t sd) { generator.seed(sd); }
    };

    /// Observes a material interactor state @param interactor_state
    template <typename interactor_state_t, typename propagator_state_t>
    DETRAY_HOST inline void operator()(state& simulator_state,
                                       interactor_state_t& interactor_state,
                                       propagator_state_t& prop_state) const {

        auto& navigation = prop_state._navigation;

        if (navigation.is_on_module()) {

            auto& stepping = prop_state._stepping;
            auto& generator = simulator_state.generator;

            const scalar_type r_theta{std::normal_distribution<scalar_type>(
                0.f, constant<scalar_type>::sqrt2 *
                         interactor_state.scattering_angle)(generator)};
            const scalar_type r_phi{std::uniform_real_distribution<scalar_type>(
                -constant<scalar_type>::pi,
                constant<scalar_type>::pi)(generator)};

            const auto dir = stepping._bound_params.dir();

            // xaxis of curvilinear plane
            const vector3 u =
                unit_vectors<vector3>().make_curvilinear_unit_u(dir);

            vector3 new_dir = axis_rotation<transform3_type>(u, r_theta)(dir);
            new_dir = axis_rotation<transform3_type>(dir, r_phi)(new_dir);

            auto& vector = stepping._bound_params.vector();

            matrix_operator().element(vector, e_bound_theta, 0u) =
                getter::theta(new_dir);
            matrix_operator().element(vector, e_bound_phi, 0u) =
                getter::phi(new_dir);

            // Flag renavigation of the current candidate
            prop_state._navigation.set_high_trust();
        }
    }
};

}  // namespace detray
