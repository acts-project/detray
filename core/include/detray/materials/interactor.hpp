/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/detail/relativistic_quantities.hpp"

namespace detray {

template <typename scalar_t>
struct interactor {

    using scalar_type = scalar_t;
    using matrix_operator = standard_matrix_operator<scalar_type>;
    using relativistic_quantities = relativistic_quantities<scalar_type>;

    template <typename material_t>
    scalar_type compute_energy_loss_bethe(
        const material_t& mat, const int pdg, const scalar_type m,
        const scalar_type qOverP, const scalar_type q = unit_constants::e) {

        const auto I = mat.material().mean_excitation_energy();
        const auto Ne = mat.material().molar_electron_density();
        // const auto thickness = slab.thickness();
        const relativistic_quantities rq(m, qOverP, q);
        const auto eps = rq.compute_epsilon(Ne, thickness);
        const auto dhalf = rq.compute_delta_half(I, Ne);
        const auto u = rq.compute_mass_term(Me);
        const auto wmax = rq.compute_WMax(m);
        // uses RPP2018 eq. 33.5 scaled from mass stopping power to linear
        // stopping power and multiplied with the material thickness to get a
        // total energy loss instead of an energy loss per length. the required
        // modification only change the prefactor which becomes identical to the
        // prefactor epsilon for the most probable value.
        const auto running = scalar_type(0.5) * std::log(u / I) +
                             scalar_type(0.5) * std::log(wmax / I) - rq.beta2 -
                             dhalf;
        return eps * running;
    }
};

}  // namespace detray