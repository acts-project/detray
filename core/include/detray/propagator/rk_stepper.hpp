/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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
#include "detray/utils/column_wise_operator.hpp"

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final
    : public base_stepper<transform3_t, constraint_t, policy_t> {

    public:
    using base_type = base_stepper<transform3_t, constraint_t, policy_t>;
    using transform3_type = transform3_t;
    using policy_type = policy_t;
    using point3 = typename transform3_type::point3;
    using vector2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    using context_type = typename magnetic_field_t::context_type;
    using matrix_operator = typename base_type::matrix_operator;
    using column_wise_op = column_wise_operator<matrix_operator>;

    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;

    DETRAY_HOST_DEVICE
    rk_stepper() {}

    struct state : public base_type::state {

        static constexpr const stepping::id id = stepping::id::e_rk;

        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type& t,
              const magnetic_field_t mag_field)
            : base_type::state(t), _magnetic_field(mag_field) {}

        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type& bound_params,
            const magnetic_field_t mag_field, const detector_t& det)
            : base_type::state(bound_params, det), _magnetic_field(mag_field) {}

        /// error tolerance
        scalar _tolerance = 1e-4;

        /// step size cutoff value
        scalar _step_size_cutoff = 1e-4;

        /// maximum trial number of RK stepping
        size_t _max_rk_step_trials = 10000;

        /// stepping data required for RKN4
        struct {
            vector3 b_first, b_middle, b_last;
            vector3 k1, k2, k3, k4;
            array_t<scalar, 4> k_qop;
        } _step_data;

        const magnetic_field_t _magnetic_field;

        // Set the local error tolerenace
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _tolerance = tol; };

        /// Update the derivative of position and direction w.r.t path length
        DETRAY_HOST_DEVICE
        inline void advance_derivative();

        /// Update the track state by Runge-Kutta-Nystrom integration.
        DETRAY_HOST_DEVICE
        inline void advance_track();

        // Update the jacobian transport from free propagation
        DETRAY_HOST_DEVICE
        inline void advance_jacobian();

        // evaulate k_n for runge kutta stepping
        DETRAY_HOST_DEVICE
        inline vector3 evaluate_k(const vector3& b_field, const int i,
                                  const scalar h, const vector3& k_prev);
    };

    /** Take a step, using an adaptive Runge-Kutta algorithm.
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation);
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"