/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/materials/interactor.hpp"

namespace detray {

template <typename interactor_t>
struct pointwise_interactor_kernel {

    using scalar_type = typename interactor_t::scalar_type;
    using vector3 = __plugin::vector3<scalar_type>;
    using interactor_type = interactor_t;
    using output_type = void;

    bool performorm_covariance_transport = true;

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

    template <typename propagator_state_t>
    pointwise_interactor_kernel(const propagator_state_t &prop_state) {
        const auto &stepping = prop_state._stepping;
        const auto &trk = stepping();

        pos = trk.pos();
        time = trk.time();
        dir = trk.dir();
        momentum = getter::norm(trk.mom());
        q = trk.charge();
        qOverP = trk.qop();
    }

    template <typename material_group_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const material_group_t &material_group, std::size_t material_index,
        const line_plane_intersection &is, propagator_state_t &prop_state,
        bool do_multiple_scattering, bool do_energy_loss) {

        const auto &mat = material_group.at(material_index);

        if (do_energy_loss) {
            e_loss = interactor_type().compute_energy_loss_bethe(
                is, mat, pdg, mass, qOverP, q);
        }

        if (performance_covariance_transport) {
            covariance_contributions(is, mat, do_multiple_scattering,
                                     do_energy_loss);
        }
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE inline convariance_contributions(
        const line_plane_intersection &is, const material_t &mat,
        bool do_multiple_scattering, bool do_energy_loss) {

        if (do_multiple_scattering) {
            // @todo: use momentum before or after energy loss in backward mode?
            scattering_angle =
                interaction_type().compute_multiple_scattering_theta0(
                    is, mat, pdg, mass, qOverP, q);
        }
        // @todo: include the radiative loss (Bremsstrahlung)
        if (do_energy_loss) {
            const scalar_type sigma_qOverP =
                interaction_type().compute_energy_loss_landau_sigma_QOverP(
                    is, mat, pdg, mass, qOverP, q);

            variance_qOverP = sigma_qOverP * sigma_qOverP;
        }
    }
};

}  // namespace detray
