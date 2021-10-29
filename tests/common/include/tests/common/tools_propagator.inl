/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "core/track.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tests/common/read_geometry.hpp"
#include "tools/line_stepper.hpp"
#include "tools/navigator.hpp"
#include "tools/propagator.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    auto [d, name_map] = read_from_csv(tml_files, host_mr);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;
    using detray_context = decltype(d)::transform_store::context;
    using detray_track = track<detray_context>;

    detray_track traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3{1., 1., 0.});
    traj.ctx = detray_context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-4;

    using detray_stepper = line_stepper<detray_track>;

    detray_stepper s;
    detray_navigator n(std::move(d));

    using detray_propagator = propagator<detray_stepper, detray_navigator>;
    detray_propagator p(std::move(s), std::move(n));

    void_track_inspector vi;

    /*auto end = */ p.propagate(traj, vi);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
