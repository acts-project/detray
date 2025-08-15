/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"

namespace detray {

template <concepts::algebra algebra_t>
struct parameter_transporter : actor {

    /// @name Type definitions for the struct
    /// @{
    using scalar_type = dscalar<algebra_t>;
    // Transformation matching this struct
    using transform3_type = dtransform3D<algebra_t>;
    // bound matrix type
    using bound_matrix_t = bound_matrix<algebra_t>;
    // free matrix type
    using free_matrix_t = free_matrix<algebra_t>;
    // Matrix type for bound to free jacobian
    using bound_to_free_matrix_t = bound_to_free_matrix<algebra_t>;
    // Matrix type for free to bound jacobian
    using free_to_bound_matrix_t = free_to_bound_matrix<algebra_t>;
    /// @}

    struct get_free_to_bound_jacobian_dlocal_dfree_kernel {
        template <typename mask_group_t, typename index_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline free_to_bound_matrix_t operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3,
            const stepper_state_t& stepping) const {
            using frame_t = typename mask_group_t::value_type::shape::
                template local_frame_type<algebra_t>;

            // Declare jacobian for free to bound coordinate transform
            free_to_bound_matrix_t jac_to_local =
                matrix::zero<free_to_bound_matrix_t>();

            detail::jacobian_engine<algebra_t>::
                template free_to_bound_jacobian_set_dlocal_dfree<frame_t>(
                    jac_to_local, trf3, stepping().pos(), stepping().dir());

            return jac_to_local;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(propagator_state_t& propagation) const {
        auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        // Geometry context for this track
        const auto& gctx = propagation._context;

        // Current Surface
        const auto sf = navigation.get_surface();

        // Bound track params of departure surface
        auto& bound_params = stepping.bound_params();

        // Covariance is transported only when the previous surface is an
        // actual tracking surface. (i.e. This disables the covariance transport
        // from curvilinear frame)
        if (!bound_params.surface_link().is_invalid()) {

            const auto full_jacobian = get_full_jacobian(propagation);

            // Calculate surface-to-surface covariance transport
            const bound_matrix_t new_cov = full_jacobian *
                                           bound_params.covariance() *
                                           matrix::transpose(full_jacobian);

            stepping.bound_params().set_covariance(new_cov);
        }

        // Convert free to bound vector
        bound_params.set_parameter_vector(
            sf.free_to_bound_vector(gctx, stepping()));

        // Set surface link
        bound_params.set_surface_link(sf.barcode());

        assert(!bound_params.is_invalid());
    }

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline free_to_bound_matrix_t get_free_to_bound_jacobian(
        propagator_state_t& propagation) const {
        const auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;
        const auto& gctx = propagation._context;
        const auto sf = navigation.get_surface();

        auto free_to_bound_jacobian = sf.template visit_mask<
            get_free_to_bound_jacobian_dlocal_dfree_kernel>(
            sf.transform(gctx), propagation._stepping);

        detail::jacobian_engine<algebra_t>::
            free_to_bound_jacobian_set_dangle_dfree_dir(free_to_bound_jacobian,
                                                        stepping().dir());

        return free_to_bound_jacobian;
    }

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline bound_matrix_t get_full_jacobian(
        propagator_state_t& propagation) const {

        const auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;

        // Geometry context for this track
        const auto& gctx = propagation._context;

        // Current Surface
        const auto sf = navigation.get_surface();

        // Bound track params of departure surface
        auto& bound_params = stepping.bound_params();

        // Previous surface
        tracking_surface prev_sf{navigation.detector(),
                                 bound_params.surface_link()};

        const bound_to_free_matrix_t bound_to_free_jacobian =
            prev_sf.bound_to_free_jacobian(gctx, bound_params);

        const bound_to_free_matrix_t free_transport_jacobian =
            stepping.transport_jacobian() * bound_to_free_jacobian;

        auto vol = navigation.get_volume();
        const auto vol_mat_ptr = vol.has_material()
                                     ? vol.material_parameters(stepping().pos())
                                     : nullptr;

        const free_matrix_t derivative =
            detail::jacobian_engine<algebra_t>::path_to_free_derivative(
                stepping().dir(), stepping.dtds(),
                stepping.dqopds(vol_mat_ptr)) *
            sf.free_to_path_derivative(gctx, stepping().pos(), stepping().dir(),
                                       stepping.dtds());

        const free_matrix_t correction_term =
            derivative + matrix::identity<free_matrix_t>();

        const free_to_bound_matrix_t free_to_bound_jacobian =
            get_free_to_bound_jacobian(propagation);

        return (free_to_bound_jacobian * correction_term) *
               free_transport_jacobian;
    }

};  // namespace detray

}  // namespace detray
