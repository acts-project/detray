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
#include "detray/propagator/composite_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/utils/curvilinear_frame.hpp"

namespace detray {

template <concepts::algebra algebra_t>
struct parameter_transporter : actor {

    /// @name Type definitions for the struct
    /// @{
    using scalar_type = dscalar<algebra_t>;
    // Transformation matching this struct
    using transform3_type = dtransform3D<algebra_t>;
    // The track parameters bound to the current sensitive/material surface
    using bound_track_parameters_type = bound_track_parameters<algebra_t>;
    // bound matrix type
    using bound_matrix_type = bound_matrix<algebra_t>;
    // Matrix type for bound to free jacobian
    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
    /// @}

    struct state {

        friend parameter_transporter;

        state() = default;

        /// Start from free track parameters
        DETRAY_HOST_DEVICE
        explicit state(const free_track_parameters<algebra_t>& free_params) {
            init(free_params);
        }

        /// Start from bound track parameters
        DETRAY_HOST_DEVICE
        explicit state(const bound_track_parameters_type& bound_params)
            : m_bound_params{bound_params} {}

        /// @returns bound track parameters - const access
        DETRAY_HOST_DEVICE
        bound_track_parameters_type& bound_params() { return m_bound_params; }

        /// Initialize the state from bound track parameters
        DETRAY_HOST_DEVICE
        void init(const bound_track_parameters_type& bound_params) {
            m_bound_params = bound_params;
        }

        /// Initialize the state from free track parameters
        DETRAY_HOST_DEVICE
        void init(const free_track_parameters<algebra_t>& free_params) {

            curvilinear_frame<algebra_t> cf(free_params);

            // Set bound track parameters
            m_bound_params.set_parameter_vector(cf.m_bound_vec);

            // A dummy covariance - should not be used
            m_bound_params.set_covariance(
                matrix::identity<bound_matrix_type>());

            // An invalid barcode - should not be used
            m_bound_params.set_surface_link(geometry::barcode{});
        }

        /// @returns bound track parameters.
        DETRAY_HOST_DEVICE
        const bound_track_parameters_type& bound_params() const {
            return m_bound_params;
        }

        /// @returns the current full Jacbian.
        DETRAY_HOST_DEVICE
        inline const bound_matrix_type& full_jacobian() const {
            return m_full_jacobian;
        }

        private:
        /// Set new full Jacbian.
        DETRAY_HOST_DEVICE
        inline void set_full_jacobian(const bound_matrix_type& jac) {
            m_full_jacobian = jac;
        }

        /// Full jacobian
        bound_matrix_type m_full_jacobian =
            matrix::identity<bound_matrix_type>();

        /// bound covariance
        bound_track_parameters_type m_bound_params{};
    };

    struct get_full_jacobian_kernel {

        template <typename mask_group_t, typename index_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline bound_matrix_type operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3,
            const bound_to_free_matrix_type& bound_to_free_jacobian,
            const material<scalar_type>* vol_mat_ptr,
            const stepper_state_t& stepping) const {

            using frame_t = typename mask_group_t::value_type::shape::
                template local_frame_type<algebra_t>;

            using jacobian_engine_t = detail::jacobian_engine<frame_t>;

            using free_matrix_t = free_matrix<algebra_t>;
            using free_to_bound_matrix_t =
                typename jacobian_engine_t::free_to_bound_matrix_type;

            // Free to bound jacobian at the destination surface
            const free_to_bound_matrix_t free_to_bound_jacobian =
                jacobian_engine_t::free_to_bound_jacobian(trf3, stepping());

            // Path correction factor
            const free_matrix_t path_correction =
                jacobian_engine_t::path_correction(
                    stepping().pos(), stepping().dir(), stepping.dtds(),
                    stepping.dqopds(vol_mat_ptr), trf3);

            const free_matrix_t correction_term =
                matrix::identity<free_matrix_t>() + path_correction;

            return free_to_bound_jacobian * correction_term *
                   stepping.transport_jacobian() * bound_to_free_jacobian;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& transporter_state,
                                       propagator_state_t& propagation) const {
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
        auto& bound_params = transporter_state.bound_params();

        // Covariance is transported only when the previous surface is an
        // actual tracking surface. (i.e. This disables the covariance transport
        // from curvilinear frame)
        if (!bound_params.surface_link().is_invalid()) {

            // Previous surface
            tracking_surface prev_sf{navigation.detector(),
                                     bound_params.surface_link()};

            const bound_to_free_matrix_type bound_to_free_jacobian =
                prev_sf.bound_to_free_jacobian(gctx, bound_params);

            auto vol = navigation.get_volume();
            const auto vol_mat_ptr =
                vol.has_material() ? vol.material_parameters(stepping().pos())
                                   : nullptr;
            transporter_state.set_full_jacobian(
                sf.template visit_mask<get_full_jacobian_kernel>(
                    sf.transform(gctx), bound_to_free_jacobian, vol_mat_ptr,
                    propagation._stepping));

            // Calculate surface-to-surface covariance transport
            const bound_matrix_type new_cov =
                transporter_state.full_jacobian() * bound_params.covariance() *
                matrix::transpose(transporter_state.full_jacobian());

            bound_params.set_covariance(new_cov);
        }

        // Convert free to bound vector
        bound_params.set_parameter_vector(
            sf.free_to_bound_vector(gctx, stepping()));

        // Set surface link
        bound_params.set_surface_link(sf.barcode());
    }
};

/// Update the stepper state after the bound track parameters were updated
template <typename algebra_t>
struct parameter_resetter : actor {

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        const parameter_transporter<algebra_t>::state& transporter_state,
        propagator_state_t& propagation) const {

        const auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        typename propagator_state_t::detector_type::geometry_context ctx{};

        // Update free params after bound params were changed by actors
        const auto sf = navigation.get_surface();
        stepping() =
            sf.bound_to_free_vector(ctx, transporter_state.bound_params());

        // Reset jacobian transport to identity matrix
        stepping.reset_transport_jacobian();
    }
};

/// Call actors that depend on the bound track parameters safely together with
/// the parameter transporter and parameter resetter
template <typename algebra_t, typename... transporter_observers>
using parameter_updater =
    composite_actor<parameter_transporter<algebra_t>, transporter_observers...,
                    parameter_resetter<algebra_t>>;

}  // namespace detray
