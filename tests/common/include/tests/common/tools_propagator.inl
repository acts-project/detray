/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/transform_store.hpp"
#include "detray/io/csv_io.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/propagator.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/read_geometry.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator) {
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

    void_track_inspector vi;

    /*auto end = */ p.propagate(traj, vi);
}
