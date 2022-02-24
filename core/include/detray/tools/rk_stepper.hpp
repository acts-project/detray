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

namespace detray {

/** Runge-Kutta-Nystrom 4th order stepper implementation */
template <typename magnetic_field_t, typename track_t,
          template <typename...> class tuple_t = dtuple,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final : public base_stepper<track_t, tuple_t> {

    public:
    using base_type = base_stepper<track_t, tuple_t>;
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using context_type = typename magnetic_field_t::context_type;

    DETRAY_HOST_DEVICE
    rk_stepper(magnetic_field_t mag_field) : _magnetic_field(mag_field) {}

    struct state : public base_type::state {
        DETRAY_HOST_DEVICE
        state(track_t& t) : base_type::state(t) {}

        // TODO: Define default step size somewhere
        scalar _step_size = 1. * unit_constants::mm;

        scalar _previous_step_size;

        // error tolerance
        scalar _tolerance = 1e-4;

        // step size cutoff value
        scalar _step_size_cutoff = 1e-4;

        // maximum trial number of RK stepping
        size_t _max_rk_step_trials = 10000;

        // accumulated path length
        scalar _path_accumulated = 0;

        // stepping data required for RKN4
        struct stepping_data {
            vector3 b_first, b_middle, b_last;
            vector3 k1, k2, k3, k4;
            array_t<scalar, 4> k_qop;
        };

        stepping_data _step_data;

        DETRAY_HOST_DEVICE
        void release_step_size() { _step_size = 1. * unit_constants::mm; }
    };

    /** Take a step, regulared by a constrained step
     *
     * @param s The state object that chaches
     * @param path_limit maximum stepsize provided by the navigator
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    DETRAY_HOST_DEVICE
    bool step(state& stepping,
              const scalar path_limit = std::numeric_limits<scalar>::max()) {
        auto& sd = stepping._step_data;

        scalar error_estimate = 0;

        // First Runge-Kutta point
        sd.b_first =
            _magnetic_field.get_field(stepping().pos(), context_type{});
        sd.k1 = evaluatek(stepping, sd.b_first, 0);

        const auto try_rk4 = [&](const scalar& h) -> bool {
            // State the square and half of the step size
            const scalar h2 = h * h;
            const scalar half_h = h * 0.5;
            auto pos = stepping().pos();
            auto dir = stepping().dir();

            // Second Runge-Kutta point
            const vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
            sd.b_middle = _magnetic_field.get_field(pos1, context_type{});
            sd.k2 = evaluatek(stepping, sd.b_middle, 1, half_h, sd.k1);

            // Third Runge-Kutta point
            sd.k3 = evaluatek(stepping, sd.b_middle, 2, half_h, sd.k2);

            // Last Runge-Kutta point
            const vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
            sd.b_last = _magnetic_field.get_field(pos2, context_type{});
            sd.k4 = evaluatek(stepping, sd.b_last, 3, h, sd.k3);

            // Compute and check the local integration error estimate
            // @Todo
            const auto err_vec = h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4);
            error_estimate =
                std::max(getter::norm(err_vec), static_cast<scalar>(1e-20));

            return (error_estimate <= stepping._tolerance);
        };

        scalar step_size_scaling = 1.;
        size_t n_step_trials = 0;

        while (!try_rk4(stepping._step_size)) {
            step_size_scaling = std::min(
                std::max(0.25,
                         std::sqrt(std::sqrt((stepping._tolerance /
                                              std::abs(2. * error_estimate))))),
                4.);

            stepping._step_size = stepping._step_size * step_size_scaling;

            // If step size becomes too small the particle remains at the
            // initial place
            if (std::abs(stepping._step_size) <
                std::abs(stepping._step_size_cutoff)) {
                // Not moving due to too low momentum needs an aborter
                printf("step size is too small. will break. \n");
                return false;
            }

            // If the parameter is off track too much or given step_size is not
            // appropriate
            if (n_step_trials > stepping._max_rk_step_trials) {
                // Too many trials, have to abort
                printf("too many rk4 trials. will break. \n");
                return false;
            }
            n_step_trials++;
        }

        stepping._step_size = std::min(stepping._step_size, path_limit);

        auto h = stepping._step_size;
        auto pos = stepping().pos();
        auto dir = stepping().dir();

        // Update the track parameters according to the equations of motion
        pos = pos + h * dir + h * h / 6. * (sd.k1 + sd.k2 + sd.k3);
        stepping().set_pos(pos);

        dir = dir + h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
        dir = vector::normalize(dir);

        stepping().set_dir(dir);

        stepping._path_accumulated += h;
        stepping._previous_step_size = h;

        // state.stepping.derivative.template head<3>() = dir;
        // state.stepping.derivative.template segment<3>(4) = sd.k4;

        return true;
    }

    DETRAY_HOST_DEVICE
    inline vector3 evaluatek(const state& stepping, const vector3& b_field,
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
    magnetic_field_t _magnetic_field;
};

}  // namespace detray