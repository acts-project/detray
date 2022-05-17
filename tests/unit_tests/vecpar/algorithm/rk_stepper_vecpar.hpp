/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#if defined(array)
#include "detray/plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#endif

#include "detray/definitions/cuda_definitions.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "vecpar/core/algorithms/parallelizable_map.hpp"
#include "vecpar/core/definitions/config.hpp"

using namespace detray;

namespace {

// type definitions
using size_type = __plugin::size_type;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using size_type = __plugin::size_type;
using mag_field_t = constant_magnetic_field<>;
using rk_stepper_t = rk_stepper<mag_field_t, free_track_parameters>;
using crk_stepper_t =
    rk_stepper<mag_field_t, free_track_parameters, constrained_step<>>;
using transform3 = __plugin::transform3<scalar>;
using matrix_operator = standard_matrix_operator<scalar>;

// geomery navigation configurations
constexpr unsigned int theta_steps = 100;
constexpr unsigned int phi_steps = 100;
constexpr unsigned int rk_steps = 100;

constexpr scalar epsilon = 1e-3;
constexpr scalar rk_tolerance = 1e-4;

// dummy navigation struct
struct nav_state {

    DETRAY_HOST_DEVICE
    scalar operator()() const { return _step_size; }

    DETRAY_HOST_DEVICE
    inline auto current_object() const -> dindex { return dindex_invalid; }

    DETRAY_HOST_DEVICE
    inline void set_full_trust() {}

    DETRAY_HOST_DEVICE
    inline void set_high_trust() {}

    DETRAY_HOST_DEVICE
    inline void set_fair_trust() {}

    DETRAY_HOST_DEVICE
    inline void set_no_trust() {}

    DETRAY_HOST_DEVICE
    inline bool abort() { return false; }

    scalar _step_size = 1. * unit_constants::mm;
};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    stepping_t _stepping;
    navigation_t _navigation;
};

}  // namespace

namespace detray {

static inline vecpar::config getConfig() {

#if !defined(__CUDA__)
    vecpar::config config{};  // let the OpenMP runtime choose
#else
    /*
     * constexpr int thread_dim = 2 * 32;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;
    vecpar::config config{block_dim, thread_dim};
     */
    vecpar::config config{1, 1};  // To match CUDA test case
#endif
    return config;
}

struct rk_stepper_algorithm
    : public vecpar::algorithm::parallelizable_mmap<
          bound_track_parameters, bound_track_parameters, vector3, transform3> {

    TARGET bound_track_parameters& map(bound_track_parameters& out_param,
                                       bound_track_parameters in_param,
                                       vector3 B, transform3 trf) override {

        mag_field_t mag_field(B);
        prop_state<crk_stepper_t::state, nav_state> propagation{
            crk_stepper_t::state(in_param, trf), nav_state{}};
        crk_stepper_t::state& crk_state = propagation._stepping;
        nav_state& n_state = propagation._navigation;

        // Decrease tolerance down to 1e-8
        crk_state.set_tolerance(rk_tolerance);

        // RK stepper and its state
        crk_stepper_t crk_stepper(mag_field);

        // Path length per turn
        scalar S = 2. * std::fabs(1. / in_param.qop()) / getter::norm(B) * M_PI;

        // Run stepper for one turn
        unsigned int max_steps = 1e4;
        for (unsigned int i = 0; i < max_steps; i++) {

            crk_state.set_constraint(S - crk_state.path_length());

            n_state._step_size = S;

            crk_stepper.step(propagation);

            if (std::abs(S - crk_state.path_length()) < 1e-6) {
                break;
            }
        }

        // Bound state after one turn propagation
        out_param = crk_stepper.bound_state(propagation, trf);
        return out_param;
    }

};  // end algorithm

}  // namespace detray