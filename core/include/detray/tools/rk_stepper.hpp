/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray tools
#include "detray/tools/base_stepper.hpp"

// detray definitions
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"

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
    using vector3 = __plugin::vector3<scalar>;
    using context_type = typename magnetic_field_t::context_type;

    DETRAY_HOST_DEVICE
    rk_stepper(const magnetic_field_t& mag_field)
        : _magnetic_field(mag_field) {}

    struct state : public base_type::state {
        DETRAY_HOST_DEVICE
        state(track_t& t) : base_type::state(t) {}

        /// error tolerance
        scalar _tolerance = 1e-4;

        /// step size cutoff value
        scalar _step_size_cutoff = 1e-4;

        /// maximum trial number of RK stepping
        size_t _max_rk_step_trials = 10000;

        /// Accumulated path length
        scalar _path_length = 0.;

        /// stepping data required for RKN4
        struct stepping_data {
            vector3 b_first, b_middle, b_last;
            vector3 k1, k2, k3, k4;
            array_t<scalar, 4> k_qop;
        };

        stepping_data _step_data;

        /// Update the track state by Runge-Kutta-Nystrom integration.
        DETRAY_HOST_DEVICE
        inline void advance_track() {
            const auto& sd = this->_step_data;
            const auto& h = this->_step_size;
            auto& track = this->_track;
            auto pos = track.pos();
            auto dir = track.dir();

            // Update the track parameters according to the equations of motion
            pos = pos + h * dir + h * h / 6. * (sd.k1 + sd.k2 + sd.k3);
            track.set_pos(pos);

            dir = dir + h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
            dir = vector::normalize(dir);
            track.set_dir(dir);

            this->_path_length += h;
        }
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
                                 navigation_state_t& navigation) {
        auto& sd = stepping._step_data;

        scalar error_estimate = 0;

        // First Runge-Kutta point
        sd.b_first =
            _magnetic_field.get_field(stepping().pos(), context_type{});
        sd.k1 = evaluate_k(stepping, sd.b_first, 0);

        const auto try_rk4 = [&](const scalar& h) -> bool {
            // State the square and half of the step size
            const scalar h2 = h * h;
            const scalar half_h = h * 0.5;
            auto pos = stepping().pos();
            auto dir = stepping().dir();

            // Second Runge-Kutta point
            const vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
            sd.b_middle = _magnetic_field.get_field(pos1, context_type{});
            sd.k2 = evaluate_k(stepping, sd.b_middle, 1, half_h, sd.k1);

            // Third Runge-Kutta point
            sd.k3 = evaluate_k(stepping, sd.b_middle, 2, half_h, sd.k2);

            // Last Runge-Kutta point
            const vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
            sd.b_last = _magnetic_field.get_field(pos2, context_type{});
            sd.k4 = evaluate_k(stepping, sd.b_last, 3, h, sd.k3);

            // Compute and check the local integration error estimate
            // @Todo
            const auto err_vec = h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4);
            error_estimate =
                std::max(getter::norm(err_vec), static_cast<scalar>(1e-20));

            return (error_estimate <= stepping._tolerance);
        };

        // Initial step size estimate
        stepping.set_step_size(navigation());
        scalar step_size_scaling = 1.;
        size_t n_step_trials = 0;

        // Adjust initial step size to integration error
        while (!try_rk4(stepping._step_size)) {

            step_size_scaling = std::min(
                std::max(0.25 * unit_constants::mm,
                         std::sqrt(std::sqrt((stepping._tolerance /
                                              std::abs(2. * error_estimate))))),
                4.);

            stepping._step_size *= step_size_scaling;

            // If step size becomes too small the particle remains at the
            // initial place
            if (std::abs(stepping._step_size) <
                std::abs(stepping._step_size_cutoff)) {
                // Not moving due to too low momentum needs an aborter
                printf("Stepper: step size is too small. will break. \n");
                // State is broken
                return navigation.abort();
            }

            // If the parameter is off track too much or given step_size is not
            // appropriate
            if (n_step_trials > stepping._max_rk_step_trials) {
                // Too many trials, have to abort
                printf("Stepper: too many rk4 trials. will break. \n");
                // State is broken
                return navigation.abort();
            }
            n_step_trials++;
        }

        // Update navigation direction
        const step::direction dir = stepping._step_size >= 0
                                        ? step::direction::e_forward
                                        : step::direction::e_backward;
        stepping.set_direction(dir);

        // Decide final step size and inform navigator
        // Not a severe change to track state expected
        if (std::abs(stepping.step_size()) <
            std::abs(
                stepping.constraints().template size<>(stepping.direction()))) {
            navigation.set_high_trust();
        }
        // Step size hit a constraint - the track state was probably changed a
        // lot
        else {
            stepping.set_step_size(
                stepping.constraints().template size<>(stepping.direction()));
            // Re-evaluate all candidates
            navigation.set_fair_trust();
        }

        // Update and check path limit
        if (not stepping.check_path_limit()) {
            printf("Stepper: Above maximal path length!\n");
            // State is broken
            return navigation.abort();
        }

        // Update track state
        stepping.advance_track();

        // state.stepping.derivative.template head<3>() = dir;
        // state.stepping.derivative.template segment<3>(4) = sd.k4;

        return true;
    }

    DETRAY_HOST_DEVICE
    inline vector3 evaluate_k(const state& stepping, const vector3& b_field,
                              const int i, const scalar h = 0.,
                              const vector3& k_prev = vector3{0, 0, 0}) {
        vector3 k_new;
        auto qop = stepping().qop();
        auto dir = stepping().dir();

        if (i == 0) {
            k_new = qop * vector::cross(dir, b_field);
        } else {
            k_new = qop * vector::cross(dir + h * k_prev, b_field);
        }

        return k_new;
    }

    private:
    const magnetic_field_t& _magnetic_field;
};

}  // namespace detray