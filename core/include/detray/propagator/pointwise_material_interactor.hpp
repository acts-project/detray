/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

template <typename interaction_t>
struct pointwise_material_interactor : actor {

    using interaction_type = interaction_t;
    using scalar_type = typename interaction_type::scalar_type;

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
        /// Evaluated varaince of qoverp
        scalar_type variance_qOverP = 0.;

        bool do_covariance_transport = true;
        bool do_energy_loss = true;
        bool do_multiple_scattering = true;
    };

    struct kernel {
        friend class pointwise_material_interactor;

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

                // Covariance update
                if (s.do_covariance_transport) {
                    if (s.do_multiple_scattering) {
                        // @todo: use momentum before or after energy loss in
                        // backward mode?
                        s.scattering_angle =
                            interaction_type()
                                .compute_multiple_scattering_theta0(
                                    is, mat, s.pdg, s.mass, qop, charge);
                    }
                    // @todo: include the radiative loss (Bremsstrahlung)
                    if (s.do_energy_loss) {
                        const scalar_type sigma_qOverP =
                            interaction_type()
                                .compute_energy_loss_landau_sigma_QOverP(
                                    is, mat, s.pdg, s.mass, qop, charge);

                        s.variance_qOverP = sigma_qOverP * sigma_qOverP;
                    }
                }
            }
            return true;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state &interactor_state, propagator_state_t &prop_state) const {

        auto &stepping = prop_state._stepping;
        auto &navigation = prop_state._navigation;

        // Do material interaction when the track is on surface
        if (navigation.is_on_module()) {

            const auto &is = *navigation.current();
            const auto &det = navigation.detector();
            const auto &surface = det->surface_by_index(is.link);
            const auto &mat_store = det->material_store();

            auto succeed = mat_store.template execute<kernel>(
                surface.material_type(), surface, is, interactor_state,
                stepping);

            if (succeed) {

                const auto &mass = interactor_state.mass;
                const scalar_type p = stepping().p();
                const auto &e_loss = interactor_state.e_loss;

                // Get new Energy
                const auto nextE =
                    std::sqrt(mass * mass + p * p) -
                    std::copysign(e_loss,
                                  static_cast<int>(navigation.direction()));

                // put particle at rest if energy loss is too large
                const auto nextP =
                    (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;

                const auto new_qop = stepping().charge() != 0.
                                         ? stepping().charge() / nextP
                                         : 1. / nextP;

                // Update bound state
                stepping._bound_params.set_qop(new_qop);
            }
        }
    }
};

}  // namespace detray