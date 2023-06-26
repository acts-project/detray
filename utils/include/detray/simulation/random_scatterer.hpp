/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/simulation/landau_distribution.hpp"
#include "detray/simulation/scattering_helper.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s).
#include <random>

namespace detray {

template <typename transform3_t>
struct random_scatterer : actor {

    using transform3_type = transform3_t;
    using matrix_operator = typename transform3_type::matrix_actor;
    using scalar_type = typename transform3_type::scalar_type;
    using vector3 = typename transform3_type::vector3;
    using interaction_type = interaction<scalar_type>;

    struct state {
        std::random_device rd{};
        std::mt19937_64 generator{rd()};

        /// The particle mass
        scalar_type mass{105.7f * unit<scalar_type>::MeV};

        /// The particle pdg
        int pdg = 13;  // default muon

        /// most probable energy loss
        scalar_type e_loss_mpv = 0.f;

        /// energy loss sigma
        scalar_type e_loss_sigma = 0.f;

        /// projected scattering angle
        scalar_type projected_scattering_angle = 0.f;

        // Simulation setup
        bool do_energy_loss = true;
        bool do_multiple_scattering = true;

        /// Constructor with seed
        ///
        /// @param sd the seed number
        state(const uint_fast64_t sd = 0u) { generator.seed(sd); }

        void set_seed(const uint_fast64_t sd) { generator.seed(sd); }
    };

    /// Material store visitor
    struct kernel {

        using scalar_type = typename interaction_type::scalar_type;
        using state = typename random_scatterer::state;

        template <typename material_group_t, typename index_t,
                  typename surface_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const material_group_t& material_group,
            const index_t& material_range,
            const intersection2D<surface_t, transform3_type>& is, state& s,
            const bound_track_parameters<transform3_type>& bound_params) const {

            const scalar qop = bound_params.qop();
            const scalar charge = bound_params.charge();

            for (const auto& mat :
                 detray::ranges::subrange(material_group, material_range)) {

                // Energy Loss
                if (s.do_energy_loss) {
                    s.e_loss_mpv =
                        interaction_type().compute_energy_loss_landau(
                            is, mat, s.pdg, s.mass, qop, charge);

                    s.e_loss_sigma =
                        interaction_type().compute_energy_loss_landau_sigma(
                            is, mat, s.pdg, s.mass, qop, charge);
                }

                // Covariance update
                if (s.do_multiple_scattering) {
                    // @todo: use momentum before or after energy loss in
                    // backward mode?
                    s.projected_scattering_angle =
                        interaction_type().compute_multiple_scattering_theta0(
                            is, mat, s.pdg, s.mass, qop, charge);
                }
            }
        }
    };

    /// Observes a material interactor state @param interactor_state
    template <typename propagator_state_t>
    DETRAY_HOST inline void operator()(state& simulator_state,
                                       propagator_state_t& prop_state) const {

        auto& navigation = prop_state._navigation;

        if (not navigation.is_on_module()) {
            return;
        }

        auto& stepping = prop_state._stepping;
        auto& bound_params = stepping._bound_params;
        const auto& is = *navigation.current();
        const auto sf = navigation.current_surface();

        sf.template visit_material<kernel>(is, simulator_state, bound_params);

        // Get the new momentum
        const auto new_mom = attenuate(
            simulator_state.e_loss_mpv, simulator_state.e_loss_sigma,
            simulator_state.mass, bound_params.p(), simulator_state.generator);

        // Update Qop
        bound_params.set_qop(bound_params.charge() / new_mom);

        // Get the new direction from random scattering
        const auto new_dir = scatter(bound_params.dir(),
                                     simulator_state.projected_scattering_angle,
                                     simulator_state.generator);

        // Update Phi and Theta
        auto& vector = stepping._bound_params.vector();
        matrix_operator().element(vector, e_bound_phi, 0u) =
            getter::phi(new_dir);
        matrix_operator().element(vector, e_bound_theta, 0u) =
            getter::theta(new_dir);

        // Flag renavigation of the current candidate
        prop_state._navigation.set_high_trust();
    }

    /// @brief Get the new momentum from the landau distribution
    template <typename generator_t>
    DETRAY_HOST inline scalar_type attenuate(const scalar_type mpv,
                                             const scalar_type sigma,
                                             const scalar_type m0,
                                             const scalar_type p0,
                                             generator_t& generator) const {

        // Get the random energy loss
        // @todo tune the scale parameters (e_loss_mpv and e_loss_sigma)
        const auto e_loss =
            landau_distribution<scalar_type>{}(generator, mpv, sigma);

        // E = sqrt(m^2 + p^2)
        const auto energy = std::sqrt(m0 * m0 + p0 * p0);
        const auto new_energy = energy - e_loss;

        auto p2 = new_energy * new_energy - m0 * m0;

        // To avoid divergence
        if (p2 < 0) {
            p2 = 1.f * unit<scalar_type>::eV * unit<scalar_type>::eV;
        }

        // p = sqrt(E^2 - m^2)
        return std::sqrt(p2);
    }

    /// @brief Scatter the direction with projected scattering angle
    ///
    /// @param dir  input direction
    /// @param projected_scattering_angle  projected scattering angle
    /// @param generator random generator
    /// @returns the new direction from random scattering
    template <typename generator_t>
    DETRAY_HOST inline vector3 scatter(
        const vector3& dir, const scalar_type projected_scattering_angle,
        generator_t& generator) const {

        // Scattering angle = sqrt(2) * projected_scattering_angle
        const auto scattering_angle =
            constant<scalar_type>::sqrt2 * projected_scattering_angle;

        return scattering_helper<transform3_type>()(dir, scattering_angle,
                                                    generator);
    }
};

}  // namespace detray
