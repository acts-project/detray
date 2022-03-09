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

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"

using namespace detray;

// type definitions
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;

using mag_field_type = constant_magnetic_field<>;
using rk_stepper_type = rk_stepper<mag_field_type, free_track_parameters>;

// geomery navigation configurations
constexpr unsigned int theta_steps = 100;
constexpr unsigned int phi_steps = 100;
constexpr unsigned int rk_steps = 100;

constexpr scalar epsilon = 1e-5;
constexpr scalar path_limit = 2 * unit_constants::m;

namespace detray {

// dummy navigation struct
struct nav_state {
    DETRAY_HOST_DEVICE scalar operator()() const {
        return 1. * unit_constants::mm;
    }
    DETRAY_HOST_DEVICE inline void set_full_trust() {}
    DETRAY_HOST_DEVICE inline void set_high_trust() {}
    DETRAY_HOST_DEVICE inline void set_fair_trust() {}
    DETRAY_HOST_DEVICE inline void set_no_trust() {}
    DETRAY_HOST_DEVICE inline bool abort() { return false; }
};

// test function for Runge-Kutta stepper
void rk_stepper_test(
    vecmem::data::vector_view<free_track_parameters>& tracks_data,
    const vector3 B);

}  // namespace detray