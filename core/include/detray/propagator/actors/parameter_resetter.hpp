/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"

namespace detray {

template <typename algebra_t>
struct parameter_resetter : actor {

    using scalar_type = dscalar<algebra_t>;

    struct state {};

    /// Mask store visitor
    struct kernel {

        // Matrix actor
        using transform3_type = dtransform3D<algebra_t>;
        using matrix_operator = dmatrix_operator<algebra_t>;

        template <typename mask_group_t, typename index_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3_type& trf3, const dindex sf_idx,
            stepper_state_t& stepping) const {

            // Note: How is it possible with "range"???
            const auto& mask = mask_group[index];

            // Reset the free vector
            stepping() = detail::bound_to_free_vector(trf3, mask,
                                                      stepping.bound_params());

            // Reset the path length
            stepping.reset_path_from_surface();

            // Reset jacobian transport to identity matrix
            stepping.reset_transport_jacobian();

            // Reset the surface index
            stepping.set_prev_sf_index(sf_idx);
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& /*resetter_state*/,
                                       propagator_state_t& propagation) const {

        const auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        using geo_cxt_t =
            typename propagator_state_t::detector_type::geometry_context;
        const geo_cxt_t ctx{};

        // Surface
        const auto sf = navigation.get_surface();

        sf.template visit_mask<kernel>(sf.transform(ctx), sf.index(), stepping);
    }
};

}  // namespace detray
