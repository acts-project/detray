/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/detail/covariance_engine.hpp"
#include "detray/propagator/navigation_policies.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename track_t,
          typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final : public base_stepper<track_t, constraint_t, policy_t> {

    public:
    using base_type = base_stepper<track_t, constraint_t, policy_t>;
    using policy_type = policy_t;
    using point3 = __plugin::point3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using context_type = typename magnetic_field_t::context_type;
    using matrix_operator = typename base_type::matrix_operator;
    using covariance_engine = typename base_type::covariance_engine;
    using vector_engine = typename covariance_engine::vector_engine;
    using transform3 = typename covariance_engine::transform3;

    DETRAY_HOST_DEVICE
    rk_stepper(const magnetic_field_t mag_field) : _magnetic_field(mag_field) {}

    struct state : public base_type::state {
        DETRAY_HOST_DEVICE
        state(const track_t& t) : base_type::state(t) {}

        DETRAY_HOST_DEVICE
        state(const bound_track_parameters& bound_params,
              const transform3& trf3)
            : base_type::state(bound_params, trf3) {}

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

        /// Set the local error tolerenace
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _tolerance = tol; };

        /// Update the derivative of position and direction w.r.t path length
        DETRAY_HOST_DEVICE
        inline void advance_derivative();

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
    };

    /// Take a step, using an adaptive Runge-Kutta algorithm.
    ///
    /// @param propagation contains the states of the rk stepper and navigator
    ///
    /// @returns the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation);

    /// Get the bound state at the surface
    ///
    /// @param propagation contains the states of the rk stepper and navigator
    /// @param trf placement of the surface
    ///
    /// @return returning the bound track parameter
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bound_track_parameters
    bound_state(propagation_state_t& propagation, const transform3& trf);

    private:
    const magnetic_field_t _magnetic_field;
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"