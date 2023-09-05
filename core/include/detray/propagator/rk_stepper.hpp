/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/navigation_policies.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename T = detray::scalar,
          template <typename> class algebra_t = ALGEBRA_PLUGIN,
          typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final
    : public base_stepper<T, algebra_t, constraint_t, policy_t> {

    public:
    using base_type = base_stepper<T, algebra_t, constraint_t, policy_t>;
    using transform3_type = dtransform3D<algebra_t<T>>;
    using point3 = dpoint3D<algebra_t<T>>;
    using vector2 = dpoint2D<algebra_t<T>>;
    using vector3 = dvector3D<algebra_t<T>>;
    using matrix_operator = typename base_type::matrix_operator;
    using mat_helper = matrix_helper<matrix_operator>;
    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;

    using policy_type = policy_t;

    DETRAY_HOST_DEVICE
    rk_stepper() {}

    struct state : public base_type::state {

        static constexpr const stepping::id id = stepping::id::e_rk;

        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type& t,
              const magnetic_field_t& mag_field)
            : base_type::state(t), _magnetic_field(mag_field) {}

        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type& bound_params,
            const magnetic_field_t& mag_field, const detector_t& det)
            : base_type::state(bound_params, det), _magnetic_field(mag_field) {}
        /// error tolerance
        scalar _tolerance{1e-4f};

        /// step size cutoff value
        scalar _step_size_cutoff{1e-4f};

        /// maximum trial number of RK stepping
        std::size_t _max_rk_step_trials{10000u};

        /// stepping data required for RKN4
        struct {
            vector3 b_first{0.f, 0.f, 0.f};
            vector3 b_middle{0.f, 0.f, 0.f};
            vector3 b_last{0.f, 0.f, 0.f};
            vector3 k1{0.f, 0.f, 0.f};
            vector3 k2{0.f, 0.f, 0.f};
            vector3 k3{0.f, 0.f, 0.f};
            vector3 k4{0.f, 0.f, 0.f};
            array_t<scalar, 4> k_qop{0.f, 0.f, 0.f, 0.f};
        } _step_data;

        /// Magnetic field view
        const magnetic_field_t _magnetic_field;

        /// Set the local error tolerenace
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _tolerance = tol; };

        /// Update the track state by Runge-Kutta-Nystrom integration.
        DETRAY_HOST_DEVICE
        inline void advance_track();

        /// Update the jacobian transport from free propagation
        DETRAY_HOST_DEVICE
        inline void advance_jacobian();

        /// evaulate k_n for runge kutta stepping
        DETRAY_HOST_DEVICE
        inline vector3 evaluate_k(const vector3& b_field, const int i,
                                  const scalar h, const vector3& k_prev);

        DETRAY_HOST_DEVICE
        inline vector3 dtds() const { return this->_step_data.k4; }
    };

    /// Take a step, using an adaptive Runge-Kutta algorithm.
    ///
    /// @return returning the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation);
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"
