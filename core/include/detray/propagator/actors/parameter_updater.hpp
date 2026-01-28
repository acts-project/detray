/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
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

namespace detray::actor {

template <concepts::algebra algebra_t>
struct parameter_transporter : base_actor {

    /// @name Type definitions for the struct
    /// @{
    using scalar_type = dscalar<algebra_t>;
    // Transformation matching this struct
    using transform3_type = dtransform3D<algebra_t>;
    // The track parameters bound to the current sensitive/material surface
    using free_track_parameters_type = free_track_parameters<algebra_t>;
    // The track parameters bound to the current sensitive/material surface
    using bound_track_parameters_type = bound_track_parameters<algebra_t>;
    // bound matrix type
    using bound_matrix_type = bound_matrix<algebra_t>;
    // free matrix type
    using free_matrix_type = free_matrix<algebra_t>;
    // Matrix type for bound to free jacobian
    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
    // Matrix type for free to bound jacobian
    using free_to_bound_matrix_type = free_to_bound_matrix<algebra_t>;
    /// @}

    struct result : public actor::result{
        /// bound track parameters of departure surface
        bound_track_parameters_type* departure_params{};
        /// bound track parameters of destination surface
        bound_track_parameters_type destination_params{};
        /// The Jacobian between the departure and destination surface
        bound_matrix_type propagation_step_jacobian{};

        /// @returns a string stream that prints the transporter result details
        DETRAY_HOST
        friend std::ostream &operator<<(std::ostream &os, const result &res) {
            os << static_cast<actor::result>(res) << std::endl;
            os << "departure params:\n" << *(res.departure_params) << std::endl;
            os << "destination params:\n" << res.destination_params << std::endl;
            return os;
        }
    }; 

    struct state {

        state() = default;

        /// Start from free track parameters
        DETRAY_HOST_DEVICE
        explicit constexpr state(const free_track_parameters_type& free_params) {
            init(free_params);
        }

        /// Start from bound track parameters
        DETRAY_HOST_DEVICE
        explicit constexpr state(const bound_track_parameters_type& bound_params) : m_bound_params{bound_params} {}

        /// Initialize the state from free track parameters
        DETRAY_HOST_DEVICE
        constexpr void init(const free_track_parameters_type& free_params) {

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
        constexpr const bound_track_parameters_type& bound_params() const {
            return m_bound_params;
        }

        /// @returns bound track parameters - const access
        DETRAY_HOST_DEVICE
        constexpr bound_track_parameters_type& bound_params() {
            return m_bound_params;
        }

        private:
        /// bound covariance
        bound_track_parameters_type m_bound_params{};
    };

    /// Filter the masks of a detector according to the local frame type
    struct select_frame {
        template <typename mask_t>
        using type = typename mask_t::local_frame;
    };

    /// Visitors to the surface local coordinate frame
    /// @{
    struct get_bound_to_free_dpos_dloc_visitor {
        template <typename frame_t>
        DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
            const frame_t& /*frame*/, const transform3_type& trf3,
            const free_track_parameters_type& params) const {

            return detail::jacobian_engine<algebra_t>::
                template bound_to_free_jacobian_submatrix_dpos_dloc<frame_t>(
                    trf3, params.pos(), params.dir());
        }
    };

    struct get_bound_to_free_dpos_dangle_visitor {
        template <typename frame_t>
        DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
            const frame_t& /*frame*/, const transform3_type& trf3,
            const free_track_parameters_type& params,
            const dmatrix<algebra_t, 3, 2>& ddir_dangle) const {

            return detail::jacobian_engine<algebra_t>::
                template bound_to_free_jacobian_submatrix_dpos_dangle<frame_t>(
                    trf3, params.pos(), params.dir(), ddir_dangle);
        }
    };

    struct get_free_to_bound_dloc_dpos_visitor {
        template <typename frame_t, typename stepper_state_t>
        DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 2, 3> operator()(
            const frame_t& /*frame*/, const transform3_type& trf3,
            const stepper_state_t& stepping) const {

            return detail::jacobian_engine<algebra_t>::
                template free_to_bound_jacobian_submatrix_dloc_dpos<frame_t>(
                    trf3, stepping().pos(), stepping().dir());
        }
    };
    /// @}

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE result operator()(state& actor_state,
                                       const propagator_state_t& propagation) const {
        const auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;

        // Result to be passed on to observing actors
        result res{actor::status::e_unknown, &actor_state.bound_params(), bound_track_parameters_type{}, bound_matrix_type{}};

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return res;
        }
        // Furthermore, there is no need to transport the initial parameters
        if (math::fabs(stepping.path_length()) == 0.f) {
            res.destination_params = actor_state.bound_params();
            res.status = actor::status::e_notify;
            return res;
        }

        // Geometry context for this track
        const auto& gctx = propagation._context;

        // Current Surface
        const auto sf = navigation.current_surface();

        // Bound track params of departure surface
        auto& bound_params = actor_state.bound_params();

        DETRAY_VERBOSE_HOST_DEVICE(
            "Actor: Transport track parameters to surface %d", sf.index());

        // Covariance is transported only when the previous surface is an
        // actual tracking surface. (i.e. This disables the covariance transport
        // from curvilinear frame).
        if (!bound_params.surface_link().is_invalid()) {

            res.propagation_step_jacobian = get_full_jacobian(propagation, bound_params);
            const bound_matrix_type& old_cov = bound_params.covariance();
            bound_matrix_type& new_cov = res.destination_params.covariance();

            detail::transport_covariance_to_bound_impl(
                old_cov, res.propagation_step_jacobian, res.destination_params.covariance());
        }

        // Convert free to bound vector
        res.destination_params.set_parameter_vector(
            sf.free_to_bound_vector(gctx, stepping()));

        // Set surface link
        res.destination_params.set_surface_link(sf.barcode());

        assert(!res.destination_params.is_invalid());

        // Notify observing actors of the results
        res.status = actor::status::e_notify;

        return res;
    }

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE constexpr bound_matrix_type get_full_jacobian(
        propagator_state_t& propagation, const bound_track_parameters_type& departure_params) const {

        // Map the surface shapes of the detector down to the common frames
        using detector_t = typename propagator_state_t::detector_type;
        using frame_registry_t =
            types::mapped_registry<typename detector_t::masks,
                                   detail::select_frame>;

        const auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;

        // Geometry context for this track
        const auto& gctx = propagation._context;

        // Our goal here is to compute the full Jacobian, which is given as:
        //
        // J_full = J_F2B * (D + I) * J_transport * J_B2F
        //
        // The transport Jacobian is given, but the free-to-bound and
        // bound-to-free matrices as well as the derivative matrix $D$ need to
        // be computed still. In order to avoid full-rank matrix
        // multiplications, we use the known substructure of the
        // aforementioned matrices in order to perform fewer operations. We
        // describe in detail the mathematics performed throughout the
        // remainder of this function.

        // Departure surface
        const tracking_surface dep_sf{navigation.detector(),
                                 departure_params.surface_link()};

        // Destination Surface
        const tracking_surface dest_sf = navigation.current_surface();

        // Free track params of departure surface
        const free_track_parameters<algebra_t> dep_free_params =
            dep_sf.bound_to_free_vector(gctx, departure_params);

        // Free track params of destination surface
        const free_track_parameters<algebra_t>& dest_free_params =
            stepping();

        // First, we will compute the bound-to-free Jacobian, which is given
        // by three sub-Jacobians. In particular, the matrix has an 8x6 shape
        // and looks like this:
        //
        //        l0  l1 phi  th q/p   t
        //  px [[  A,  A,  B,  B,  0,  0],
        //  py  [  A,  A,  B,  B,  0,  0],
        //  pz  [  A,  A,  B,  B,  0,  0],
        //   t  [  0,  0,  0,  0,  0,  1],
        //  dx  [  0,  0,  C,  C,  0,  0],
        //  dy  [  0,  0,  C,  C,  0,  0],
        //  dz  [  0,  0,  0,  C,  0,  0],
        // q/p  [  0,  0,  0,  0,  1,  0]]
        //
        // Note that we thus only have 17 out of 48 non-trivial matrix
        // elements.
        //
        // In this matrix, A represents the d(pos)/d(loc) submatrix, B is the
        // d(pos)/d(angle) submatrix, and C is the d(dir)/d(angle) submatrix,
        // all of which are 3x2 in size. Also, A and B depend on the frame
        // type while C is computed in the same way for all frames. Finally,
        // note that submatrix B is the zero matrix for most frame types.
        const auto& dep_trf3 = dep_sf.transform(gctx);
        const dmatrix<algebra_t, 3, 2> b2f_dpos_dloc =
            types::visit<frame_registry_t, get_bound_to_free_dpos_dloc_visitor>(
                dep_sf.shape_id(), dep_trf3, dep_free_params);

        const dmatrix<algebra_t, 3, 2> b2f_ddir_dangle =
            detail::jacobian_engine<algebra_t>::
                bound_to_free_jacobian_submatrix_ddir_dangle(departure_params);

        const dmatrix<algebra_t, 3, 2> b2f_dpos_dangle =
            types::visit<frame_registry_t,
                         get_bound_to_free_dpos_dangle_visitor>(
                dep_sf.shape_id(), dep_trf3, dep_free_params, b2f_ddir_dangle);

        // Next, we compute the derivative which is defined as the outer
        // product of the two 8x1 vectors representing the path to free
        // derivative and the free to path derivative. Thus, it is the product
        // of the following two vectors:
        //
        //  px [[ A],  [[ B],^T
        //  py  [ A],   [ B],
        //  pz  [ A],   [ B],
        //   t  [ 0]]   [ 0],
        //  dx  [ A],   [ C],
        //  dy  [ A],   [ C],
        //  dz  [ A],   [ C],
        // q/p  [ A],   [ 0]]
        //
        // Where A is frame independent and non-zero, B is frame-dependent and
        // non-zero, while C is frame-dependent and non-zero only for line
        // frames.
        const dpoint3D<algebra_t> dest_glob_pos{dest_free_params.pos()};
        const dvector3D<algebra_t> dest_glob_dir{dest_free_params.dir()};

        auto vol = navigation.current_volume();
        const auto vol_mat_ptr = vol.has_material()
                                     ? vol.material_parameters(dest_glob_pos)
                                     : nullptr;

        const auto path_to_free_derivative =
            detail::jacobian_engine<algebra_t>::path_to_free_derivative(
                dest_glob_dir, stepping.dtds(),
                stepping.dqopds(vol_mat_ptr));

        const auto free_to_path_derivative = dest_sf.free_to_path_derivative(
            gctx, dest_glob_pos, dest_glob_dir, stepping.dtds());

        // Now, we compute the free-to-bound Jacobian which is of size 6x8 and
        // has the following structure:
        //
        //        px  py  pz   t  dx  dy  dz q/p
        //  l0 [[  A,  A,  A,  0,  0,  0,  0,  0],
        //  l1  [  A,  A,  A,  0,  0,  0,  0,  0],
        // phi  [  0,  0,  0,  0,  B,  B,  0,  0],
        //  th  [  0,  0,  0,  0,  B,  B,  B,  0],
        // q/p  [  0,  0,  0,  0,  0,  0,  0,  1],
        //   t  [  0,  0,  0,  1,  0,  0,  0,  0]]
        //
        // Thus, the number of non-trivial elements is only 11 out of the 48
        // matrix elements. Also, submatrix A depends on the frame type while
        // submatrix B is the same for all frame types.
        const dmatrix<algebra_t, 2, 3> f2b_dloc_dpos =
            types::visit<frame_registry_t, get_free_to_bound_dloc_dpos_visitor>(
                dest_sf.shape_id(), dest_sf.transform(gctx), stepping);

        const dmatrix<algebra_t, 2, 3> f2b_dangle_ddir =
            detail::jacobian_engine<algebra_t>::
                free_to_bound_jacobian_submatrix_dangle_ddir(dest_glob_dir);

        // Finally, we can use our Sympy-generated full Jacobian computation
        // and return its result.
        bound_matrix_type full_jacobian;

        detail::update_full_jacobian_impl(
            stepping.transport_jacobian(), b2f_dpos_dloc, b2f_ddir_dangle,
            b2f_dpos_dangle, path_to_free_derivative, free_to_path_derivative,
            f2b_dloc_dpos, f2b_dangle_ddir, full_jacobian);

        return full_jacobian;
    }
};

template <concepts::algebra algebra_t>
struct parameter_setter : base_actor {

    // The track parameters bound to the current sensitive/material surface
    using bound_track_parameters_type = bound_track_parameters<algebra_t>;
    // bound matrix type
    using bound_matrix_type = bound_matrix<algebra_t>;

    struct noise_cfg {
        /// Percentage of total track path to assume as accumulated error
        float accumulated_error{0.001f};
        /// Number of standard deviations to assume to model the scattering
        /// noise
        int n_stddev{2};
        /// Estimate mask tolerance for navigation to for compensate scattering
        bool estimate_scattering_noise{true};
    };

    struct state {

        friend parameter_setter;

        /// Default construction
        state() = default;

        /// Build from propagation configuration set by the user
        DETRAY_HOST_DEVICE
        explicit constexpr state(const propagation::config& cfg, bound_matrix_type* full_jac = nullptr)
            :  m_full_jacobian(full_jac), 
              m_cfg{cfg.navigation.accumulated_error, cfg.navigation.n_scattering_stddev, cfg.navigation.estimate_scattering_noise} {}

        /// @return access to the noise estimation configuration
        const noise_cfg& noise_estimation_cfg() const {
            return m_cfg;
        }

        /// @returns true if the full Jacobian matrix should be assembled.
        DETRAY_HOST_DEVICE
        constexpr bool has_full_jacobian() const {
            return m_full_jacobian != nullptr;
        }

        /// @returns the current full Jacbian.
        DETRAY_HOST_DEVICE
        constexpr const bound_matrix_type& full_jacobian() const {
            return *m_full_jacobian;
        }

        private:
        /// Set new full Jacbian.
        DETRAY_HOST_DEVICE
        constexpr void set_full_jacobian(const bound_matrix_type& jac) {
            *m_full_jacobian = jac;
        }

        /// Full jacobian for up to the current destination surface
        bound_matrix_type* m_full_jacobian{nullptr};
        /// Configuration for the niose estimation
        noise_cfg m_cfg{};
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& setter_state, 
                const typename parameter_transporter<algebra_t>::result& res,
                                       propagator_state_t& propagation) const {

        using scalar_t = dscalar<algebra_t>;

        // Update the track parameters if necessary
        if (res.status != actor::status::e_success) {
            return;
        }

        const auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        assert(navigation.is_on_surface());

        // If the navigation trust level has been reduced, the track parameters
        // might have been updated. Otherwise do nothing before prop. starts
        DETRAY_VERBOSE_HOST_DEVICE("Actor: Update the track parameters");

        // Update the bound parameters of the propagation flow
        bound_track_parameters_type& bound_params = *(res.departure_params);
        bound_params = res.destination_params;

        // No need to update the free params on the seed: only the cov. changed
        if (math::fabs(stepping.path_length()) > 0.f) {
            // Update free params after bound params were changed by actors
            const tracking_surface dest_sf{navigation.detector(), bound_params.surface_link()};
            stepping() =
                dest_sf.bound_to_free_vector(propagation._context, bound_params);
            assert(!stepping().is_invalid());

            // Update the full Jacobian, if required
            if (setter_state.has_full_jacobian()) {
                const auto aggregate_full_jacobian =
                    res.propagation_step_jacobian * setter_state.full_jacobian();
                setter_state.set_full_jacobian(aggregate_full_jacobian);
            }

            // Reset transport Jacobian to identity matrix
            stepping.reset_transport_jacobian();
        }

        // Track pos/dir is not known precisely: adjust navigation tolerances
        const noise_cfg& cfg = setter_state.noise_estimation_cfg();
        if (cfg.estimate_scattering_noise) {
            detail::estimate_external_mask_tolerance(
                bound_params, propagation,
                static_cast<scalar_t>(cfg.n_stddev),
                cfg.accumulated_error);
        }
    }
};

/// Call actors that depend on the bound track parameters safely together with
/// the parameter transporter and parameter setter
template <typename algebra_t, typename... transporter_observers>
using parameter_updater =
    composite_actor<parameter_transporter<algebra_t>, transporter_observers...,
                    parameter_setter<algebra_t>>;

}  // namespace detray::actor