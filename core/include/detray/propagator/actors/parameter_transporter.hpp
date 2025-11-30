/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "./codegen/covariance_transport.hpp"
#include "./codegen/full_jacobian.hpp"
#include "algebra/utils/approximately_equal.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/utils/log.hpp"
#include "detray/utils/type_registry.hpp"

namespace detray {

namespace detail {

/// Filter the masks of a detector according to the local frame type
struct select_frame {
    template <typename mask_t>
    using type = typename mask_t::local_frame;
};

}  // namespace detail

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

    struct get_bound_to_free_dpos_dloc_visitor {
        template <typename frame_t>
        DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
            const frame_t& /*frame*/, const transform3_type& trf3,
            const free_track_parameters<algebra_t>& params) const {

            return detail::jacobian_engine<algebra_t>::
                template bound_to_free_jacobian_submatrix_dpos_dloc<frame_t>(
                    trf3, params.pos(), params.dir());
        }
    };

    struct get_bound_to_free_dpos_dangle_visitor {
        template <typename frame_t>
        DETRAY_HOST_DEVICE constexpr dmatrix<algebra_t, 3, 2> operator()(
            const frame_t& /*frame*/, const transform3_type& trf3,
            const free_track_parameters<algebra_t>& params,
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

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(propagator_state_t& propagation) const {
        auto& stepping = propagation._stepping;
        const auto& navigation = propagation._navigation;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }
        DETRAY_VERBOSE_HOST_DEVICE(
            "Transport track parameters to current surface");

        // Geometry context for this track
        const auto& gctx = propagation._context;

        // Current Surface
        const auto sf = navigation.current_surface();

        // Bound track params of departure surface
        auto& bound_params = stepping.bound_params();

        // Covariance is transported only when the previous surface is an
        // actual tracking surface. (i.e. This disables the covariance transport
        // from curvilinear frame)
        if (!bound_params.surface_link().is_invalid()) {
            const auto full_jacobian = get_full_jacobian(propagation);
            const bound_matrix_t old_cov = stepping.bound_params().covariance();

            detail::transport_covariance_to_bound_impl(
                old_cov, full_jacobian, stepping.bound_params().covariance());
        }

        // Convert free to bound vector
        bound_params.set_parameter_vector(
            sf.free_to_bound_vector(gctx, stepping()));

        // Set surface link
        bound_params.set_surface_link(sf.barcode());

        assert(!bound_params.is_invalid());
    }

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE constexpr bound_matrix_t get_full_jacobian(
        propagator_state_t& propagation) const {

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

        // Current Surface
        const auto sf = navigation.current_surface();

        // Bound track params of departure surface
        auto& bound_params = stepping.bound_params();

        // Previous surface
        tracking_surface prev_sf{navigation.detector(),
                                 bound_params.surface_link()};

        // Free track params of departure surface
        const free_track_parameters<algebra_t> free_params =
            prev_sf.bound_to_free_vector(gctx, bound_params);

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
        const auto& prev_trf3 = prev_sf.transform(gctx);
        const dmatrix<algebra_t, 3, 2> b2f_dpos_dloc =
            types::visit<frame_registry_t, get_bound_to_free_dpos_dloc_visitor>(
                prev_sf.shape_id(), prev_trf3, free_params);

        const dmatrix<algebra_t, 3, 2> b2f_ddir_dangle =
            detail::jacobian_engine<algebra_t>::
                bound_to_free_jacobian_submatrix_ddir_dangle(bound_params);

        const dmatrix<algebra_t, 3, 2> b2f_dpos_dangle =
            types::visit<frame_registry_t,
                         get_bound_to_free_dpos_dangle_visitor>(
                prev_sf.shape_id(), prev_trf3, free_params, b2f_ddir_dangle);

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
        auto vol = navigation.current_volume();
        const auto vol_mat_ptr = vol.has_material()
                                     ? vol.material_parameters(stepping().pos())
                                     : nullptr;

        const auto path_to_free_derivative =
            detail::jacobian_engine<algebra_t>::path_to_free_derivative(
                stepping().dir(), stepping.dtds(),
                stepping.dqopds(vol_mat_ptr));

        const auto free_to_path_derivative = sf.free_to_path_derivative(
            gctx, stepping().pos(), stepping().dir(), stepping.dtds());

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
                sf.shape_id(), sf.transform(gctx), propagation._stepping);

        const dmatrix<algebra_t, 2, 3> f2b_dangle_ddir =
            detail::jacobian_engine<algebra_t>::
                free_to_bound_jacobian_submatrix_dangle_ddir(stepping().dir());

        // Finally, we can use our Sympy-generated full Jacobian computation
        // and return its result.
        bound_matrix_t full_jacobian;

        detail::update_full_jacobian_impl(
            stepping.transport_jacobian(), b2f_dpos_dloc, b2f_ddir_dangle,
            b2f_dpos_dangle, path_to_free_derivative, free_to_path_derivative,
            f2b_dloc_dpos, f2b_dangle_ddir, full_jacobian);

        return full_jacobian;
    }
};  // namespace detray

}  // namespace detray
