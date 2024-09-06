/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/policies.hpp"
#include "detray/propagator/base_stepper.hpp"

namespace detray {

/// Straight line stepper implementation
template <typename algebra_t, typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy,
          typename inspector_t = stepping::void_inspector>
class line_stepper final
    : public base_stepper<algebra_t, constraint_t, policy_t, inspector_t> {

    public:
    using base_type =
        base_stepper<algebra_t, constraint_t, policy_t, inspector_t>;

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;
    using matrix_operator = typename base_type::matrix_operator;
    template <std::size_t ROWS, std::size_t COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

    struct state : public base_type::state {
        static constexpr const stepping::id id = stepping::id::e_linear;

        using base_state = typename base_type::state;

        /// Import base class constructors
        using base_state::base_state;

        /// Update the track state in a straight line.
        DETRAY_HOST_DEVICE
        inline void advance_track() {
            auto& track = this->_track;
            track.set_pos(track.pos() + track.dir() * this->_step_size);

            this->_path_length += this->_step_size;
            this->_s += this->_step_size;
        }

        DETRAY_HOST_DEVICE
        inline void advance_jacobian() {

            // The step transport matrix in global coordinates
            free_matrix<algebra_t> D =
                matrix_operator().template identity<e_free_size, e_free_size>();

            // d(x,y,z)/d(n_x,n_y,n_z)
            matrix_type<3, 3> dxdn =
                this->_step_size * matrix_operator().template identity<3, 3>();
            matrix_operator().template set_block<3, 3>(D, dxdn, e_free_pos0,
                                                       e_free_dir0);

            /// NOTE: Let's skip the element for d(time)/d(qoverp) for the
            /// moment..

            this->_jac_transport = D * this->_jac_transport;
        }

        DETRAY_HOST_DEVICE
        inline vector3_type dtds() const { return {0.f, 0.f, 0.f}; }

        DETRAY_HOST_DEVICE
        inline scalar_type dqopds() const { return 0.f; }
    };

    /// Take a step, regulared by a constrained step
    ///
    /// @param stepping The state object of a stepper
    /// @param navigation The state object of a navigator
    /// @param max_step_size Maximal distance for this step
    ///
    /// @return returning the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation,
                                 const stepping::config& cfg = {}) const {
        // Get stepper and navigator states
        state& stepping = propagation._stepping;
        auto& navigation = propagation._navigation;

        if (stepping._step_size == 0.f) {
            stepping._step_size = navigation();
        } else if (stepping._step_size > 0) {
            stepping._step_size = math::min(stepping._step_size, navigation());
        } else {
            stepping._step_size = math::max(stepping._step_size, navigation());
        }

        // Escape the initialized state
        stepping._initialized = false;

        // Update navigation direction
        const step::direction step_dir = stepping._step_size >= 0.f
                                             ? step::direction::e_forward
                                             : step::direction::e_backward;
        stepping.set_direction(step_dir);

        // Check constraints
        if (math::fabs(stepping.step_size()) >
            math::fabs(
                stepping.constraints().template size<>(stepping.direction()))) {
            // Run inspection before step size is cut
            stepping.run_inspector(cfg, "Before constraint: ");

            stepping.set_step_size(
                stepping.constraints().template size<>(stepping.direction()));
        }

        // Update track state
        stepping.advance_track();

        // Advance jacobian transport
        stepping.advance_jacobian();

        // Call navigation update policy
        typename line_stepper::policy_type{}(stepping.policy_state(),
                                             propagation);

        // Run inspection if needed
        stepping.run_inspector(cfg, "Step complete: ");

        return true;
    }
};

}  // namespace detray
