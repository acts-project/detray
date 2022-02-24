/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/transform_store.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/io/csv_io.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/propagator.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/helix_gun.hpp"
#include "tests/common/tools/read_geometry.hpp"

constexpr scalar epsilon = 1e-4;

// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator_line_stepper) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    // auto [d, name_map] = read_from_csv(tml_files, host_mr);
    auto d = create_toy_geometry(host_mr);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;
    using detray_track = free_track_parameters;

    __plugin::point3<scalar> pos{0., 0., 0.};
    __plugin::vector3<scalar> mom{1., 1., 0.};
    detray_track traj(pos, 0, mom, -1);

    using detray_stepper = line_stepper<detray_track>;

    detray_stepper s;
    detray_navigator n(std::move(d));

    using detray_propagator = propagator<detray_stepper, detray_navigator>;
    detray_propagator p(std::move(s), std::move(n));

    void_propagator_inspector vi;

    /*auto end = */ p.propagate(traj, vi);
}

struct helix_inspector {

    helix_inspector(const helix_gun& helix) : _helix(helix) {}

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t& /*navigation*/,
                                       const stepper_state_t& stepping) {
        auto pos = stepping().pos();
        auto true_pos = _helix(stepping._path_accumulated);

        auto relative_error = 1 / stepping._path_accumulated * (pos - true_pos);

        EXPECT_NEAR(getter::norm(relative_error), 0, epsilon);
    }

    helix_gun _helix;
};

struct print_inspector {

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t& navigation,
                                       const stepper_state_t& stepping) {
        std::stringstream stream;

        stream << std::left << std::setw(30);
        switch (static_cast<int>(navigation.status())) {
            case -3:
                stream << "status: on_target";
                break;
            case -2:
                stream << "status: abort";
                break;
            case -1:
                stream << "status: unknowm";
                break;
            case 0:
                stream << "status: towards_surface";
                break;
            case 1:
                stream << "status: on_surface";
                break;
        };

        if (navigation.volume() == dindex_invalid) {
            stream << "volume: " << std::setw(10) << "invalid";
        } else {
            stream << "volume: " << std::setw(10) << navigation.volume();
        }

        if (navigation.on_object() == dindex_invalid) {
            stream << "surface: " << std::setw(14) << "invalid";
        } else {
            stream << "surface: " << std::setw(14) << navigation.on_object();
        }

        stream << "step_size: " << std::setw(10)
               << stepping._previous_step_size;

        std::cout << stream.str() << std::endl;
    }
};

struct combined_inspector {
    helix_inspector _hi;
    print_inspector _pi;

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t& navigation,
                                       const stepper_state_t& stepping) {
        _hi(navigation, stepping);
        _pi(navigation, stepping);
    }
};

TEST(ALGEBRA_PLUGIN, propagator_rk_stepper) {

    using namespace detray;
    using namespace __plugin;
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;
    using detray_b_field = constant_magnetic_field<>;
    using detray_stepper = rk_stepper<detray_b_field, free_track_parameters>;
    using detray_propagator = propagator<detray_stepper, detray_navigator>;

    // Constant magnetic field
    vector3 B{0, 0, 2 * unit_constants::T};
    detray_b_field b_field(B);

    // Set initial position and momentum of tracks
    const point3 pos{0., 0., 0.};
    const vector3 mom{1., 1., 0.};

    detray_stepper s(b_field);
    detray_navigator n(std::move(d));
    detray_propagator p(std::move(s), std::move(n));

    free_track_parameters traj(pos, 0, mom, -1);
    helix_gun helix(traj, B);

    combined_inspector ci{helix_inspector(helix), print_inspector{}};

    p.propagate(traj, ci);
}
