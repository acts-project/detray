/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/axis_rotation.hpp"

namespace detray {

template <typename transform3_t>
struct pointwise_material_interactor : actor {
    using transform3_type = transform3_t;
    using matrix_operator = typename transform3_t::matrix_actor;
    using scalar_type = typename matrix_operator::scalar_type;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using interaction_type = interaction<scalar_type>;
    using vector3 = typename transform3_t::vector3;

    struct state {
        using vector3 = __plugin::vector3<scalar>;

        /// The particle mass
        scalar_type mass = 105.7 * unit_constants::MeV;
        /// The particle pdg
        int pdg = 13;  // default muon

        /// Evaluated energy loss
        scalar_type e_loss = 0.;
        /// Evaluated multiple scattering angle
        scalar_type scattering_angle = 0.;
        /// Evaluated sigma of qoverp
        scalar_type sigma_qop = 0.;

        bool do_covariance_transport = true;
        bool do_energy_loss = true;
        bool do_multiple_scattering = true;

        DETRAY_HOST_DEVICE
        void reset() {
            e_loss = 0;
            scattering_angle = 0;
            sigma_qop = 0;
        }
    };

    struct kernel {

        using output_type = bool;
        using scalar_type = typename interaction_type::scalar_type;
        using state = typename pointwise_material_interactor::state;

        template <typename material_group_t, typename surface_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const material_group_t &material_group, const surface_t &surface,
            const line_plane_intersection &is, state &s,
            const stepper_state_t &stepping) const {

            const auto &material_range = surface.material_range();

            const scalar qop = stepping().qop();
            const scalar charge = stepping().charge();

            for (const auto &mat : range(material_group, material_range)) {

                // Energy Loss
                if (s.do_energy_loss) {
                    s.e_loss = interaction_type().compute_energy_loss_bethe(
                        is, mat, s.pdg, s.mass, qop, charge);
                }

                // @todo: include the radiative loss (Bremsstrahlung)
                if (s.do_energy_loss && s.do_covariance_transport) {
                    s.sigma_qop = interaction_type()
                                      .compute_energy_loss_landau_sigma_QOverP(
                                          is, mat, s.pdg, s.mass, qop, charge);
                }

                // Covariance update
                if (s.do_multiple_scattering) {
                    // @todo: use momentum before or after energy loss in
                    // backward mode?
                    s.scattering_angle =
                        interaction_type().compute_multiple_scattering_theta0(
                            is, mat, s.pdg, s.mass, qop, charge);
                }
            }
            return true;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state &interactor_state, propagator_state_t &prop_state) const {

        interactor_state.reset();

        auto &stepping = prop_state._stepping;
        auto &navigation = prop_state._navigation;

        // Do material interaction when the track is on surface
        if (navigation.is_on_module()) {

            const auto &is = *navigation.current();
            const auto &det = navigation.detector();
            const auto &surface = det->surface_by_index(is.index);
            const auto &mat_store = det->material_store();

            auto succeed = mat_store.template execute<kernel>(
                surface.material_type(), surface, is, interactor_state,
                stepping);

            if (succeed) {

                update_qop(stepping, interactor_state,
                           static_cast<int>(navigation.direction()));

                update_angle(stepping, interactor_state,
                             static_cast<int>(navigation.direction()));
            }
        }
    }

    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline void update_qop(stepper_state_t &stepping,
                                              const state &s, int sign) const {

        // Update bound qop value
        if (s.do_energy_loss) {
            const auto &m = s.mass;
            const scalar_type p = stepping().p();

            // Get new Energy
            const auto nextE =
                std::sqrt(m * m + p * p) - std::copysign(s.e_loss, sign);

            // put particle at rest if energy loss is too large
            const auto nextP =
                (m < nextE) ? std::sqrt(nextE * nextE - m * m) : 0;

            const auto new_qop = stepping().charge() != 0.
                                     ? stepping().charge() / nextP
                                     : 1. / nextP;

            stepping._bound_params.set_qop(new_qop);

            // Update bound qop varaince
            if (s.do_covariance_transport) {
                auto &covariance = stepping._bound_params.covariance();

                const scalar_type variance_qop = s.sigma_qop * s.sigma_qop;

                matrix_operator().element(covariance, e_bound_qoverp,
                                          e_bound_qoverp) +=
                    std::copysign(variance_qop, sign);
            }
        }
    }

    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline void update_angle(stepper_state_t &stepping,
                                                const state &s,
                                                int sign) const {
        if (s.do_covariance_transport) {

            auto &covariance = stepping._bound_params.covariance();

            // variance of scattering angle
            const scalar_type var_scattering_angle =
                s.scattering_angle * s.scattering_angle;

            const auto dir = stepping._bound_params.dir();

            matrix_operator().element(covariance, e_bound_phi, e_bound_phi) +=
                sign * var_scattering_angle /
                (dir[0] * dir[0] + dir[1] * dir[1]);

            matrix_operator().element(covariance, e_bound_theta,
                                      e_bound_theta) +=
                sign * var_scattering_angle;
        }
    }
};

}  // namespace detray