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
    using matrix_operator = standard_matrix_operator<scalar>;

    struct state {
        using scalar_type = typename interaction_type::scalar_type;
        using vector3 = __plugin::vector3<scalar>;

        template <typename propagator_state_t>
        DETRAY_HOST_DEVICE void set_state(propagator_state_t prop_state) {
            const auto &stepping = prop_state._stepping;
            const auto &trk = stepping();

            pos = trk.pos();
            time = trk.time();
            dir = trk.dir();
            momentum = getter::norm(trk.mom());
            q = trk.charge();
            qOverP = trk.qop();
        }

        DETRAY_HOST_DEVICE scalar_type energy() {
            return std::sqrt(mass * mass + momentum * momentum);
        }

        /// The particle position at the interaction.
        vector3 pos = {0., 0., 0};
        /// The particle time at the interaction.
        scalar_type time = 0.0;
        /// The particle direction at the interaction.
        vector3 dir = {0., 0., 0};
        /// The particle momentum at the interaction
        scalar_type momentum;
        /// The particle charge
        scalar_type q;
        /// The particle q/p at the interaction
        scalar_type qOverP;
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

        template <typename material_group_t, typename surface_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const material_group_t &material_group, const surface_t &surface,
            const line_plane_intersection &is, state &s) const {

            const auto &material_range = surface.material_range();

            for (const auto &mat : range(material_group, material_range)) {

                // Energy Loss
                if (s.do_energy_loss) {
                    s.e_loss = interaction_type().compute_energy_loss_bethe(
                        is, mat, s.pdg, s.mass, s.qOverP, s.q);
                }

                // Covariance update
                if (s.do_covariance_transport) {
                    if (s.do_multiple_scattering) {
                        // @todo: use momentum before or after energy loss in
                        // backward mode?
                        s.scattering_angle =
                            interaction_type()
                                .compute_multiple_scattering_theta0(
                                    is, mat, s.pdg, s.mass, s.qOverP, s.q);
                    }
                    // @todo: include the radiative loss (Bremsstrahlung)
                    if (s.do_energy_loss) {
                        const scalar_type sigma_qOverP =
                            interaction_type()
                                .compute_energy_loss_landau_sigma_QOverP(
                                    is, mat, s.pdg, s.mass, s.qOverP, s.q);

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

        interactor_state.set_state(prop_state);

        auto &stepping = prop_state._stepping;
        auto &navigation = prop_state._navigation;

        // Do material interaction when the track is on surface
        if (navigation.is_on_module()) {

            const auto &is = *navigation.current();
            const auto &det = navigation.detector();
            const auto &surface = det->surface_by_index(is.link);
            const auto &mat_store = det->material_store();

            auto succeed = mat_store.template execute<kernel>(
                surface.material_type(), surface, is, interactor_state);

            if (succeed) {

                const auto &mass = interactor_state.mass;
                const auto &mom = interactor_state.momentum;
                const auto &e_loss = interactor_state.e_loss;

                // Get new Energy
                const auto nextE =
                    std::sqrt(mass * mass + mom * mom) -
                    std::copysign(e_loss,
                                  static_cast<int>(navigation.direction()));

                // put particle at rest if energy loss is too large
                const auto nextP =
                    (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;

                const auto new_qop = interactor_state.q != 0.
                                         ? interactor_state.q / nextP
                                         : 1. / nextP;

                // Update momentum
                stepping().set_qop(new_qop);  // Skip free track params update?

                auto &bvec = stepping._bound_params.vector();
                // auto& bcov = stepping._bound_params.covariance();
                matrix_operator().element(bvec, e_bound_qoverp, 0) = new_qop;
            }
        }
    }
};

}  // namespace detray