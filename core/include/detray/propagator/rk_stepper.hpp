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

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename track_t,
          typename constraint_t = unconstrained_step,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final : public base_stepper<track_t, constraint_t> {

    public:
    using base_type = base_stepper<track_t, constraint_t>;
    using point3 = __plugin::point3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;
    using context_type = typename magnetic_field_t::context_type;

    DETRAY_HOST_DEVICE
    rk_stepper(const magnetic_field_t mag_field) : _magnetic_field(mag_field) {}

    struct state : public base_type::state {
        DETRAY_HOST_DEVICE
        state(track_t& t) : base_type::state(t) {}

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

        /// Update the track state by Runge-Kutta-Nystrom integration.
        DETRAY_HOST_DEVICE
        inline void advance_track();

        // Update the jacobian transport from free propagation
        DETRAY_HOST_DEVICE
        inline void advance_jacobian();

        DETRAY_HOST_DEVICE
        inline vector3 evaluate_k(const vector3& b_field, const int i,
                                  const scalar h, const vector3& k_prev);
    };

    /** Take a step, using an adaptive Runge-Kutta algorithm.
     *
     * @param stepping The state object of a stepper
     * @param navigation The state object of a navigator
     * @param max_step_size Maximal distance for this step
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    template <typename navigation_state_t>
    DETRAY_HOST_DEVICE bool step(state& stepping,
                                 navigation_state_t& navigation);

    private:
    const magnetic_field_t _magnetic_field;
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"