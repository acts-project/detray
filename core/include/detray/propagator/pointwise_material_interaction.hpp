/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/materials/interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/navigator.hpp"

namespace detray {

template <scalar_t>
struct pointwise_material_interactor : actor {

    using scalar_type = scalar_t;
    using vector3 = __plugin::vector3<scalar_type>;
    using interactor_type = interactor<scalar_type>;

    struct state {

        /// The surface id
        // const std::size_t surface_id;

        /// The particle position at the interaction.
        // const vector3 pos = {0., 0., 0};
        /// The particle time at the interaction.
        // const scalar_type time = 0.0;
        /// The particle direction at the interaction.
        // const vector3 dir = {0., 0., 0};
        /// The particle momentum at the interaction
        // const scalar momentum;
        /// The particle charge
        // const scalar q;
        /// The particle q/p at the interaction
        // const scalar qOverP;
        /// The particle mass
        // const scalar mass;
        /// The particle pdg
        // const int pdg;

        /// The covariance transport decision at the interaction
        const bool perform_covariance_transport = true;
        /// The navigation direction
        const navigation::direction nav;

        /// The effective, passed material properties including the path
        /// correction.
        /// MaterialSlab slab;
        /// The path correction factor due to non-zero incidence on the surface.
        // scalar_type pathCorrection = 0.;
        /// Expected phi variance due to the interactions.
        // scalar_type variancePhi = 0.;
        /// Expected theta variance due to the interactions.
        // scalar_type varianceTheta = 0.;
        /// Expected q/p variance due to the interactions.
        // scalar_type varianceQoverP = 0.;
        /// The energy change due to the interaction.
        scalar_type e_loss = 0.;
        /// The momentum after the interaction
        // scalar_type nextP;
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &interactor_state,
                                       propagator_state_t &prop_state,
                                       bool energy_loss,
                                       bool multiple_scattering) const {
        const auto is = *prop_state._navigation.current();

        /// Evaluate interaction
        if (energy_loss) {
            // e_loss = ;
        }
    }
};

}  // namespace detray