/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/navigation/direct_navigator.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/utils/detectors/build_toy_detector.hpp"
#include "detray/test/utils/detectors/build_wire_chamber.hpp"
#include "detray/test/utils/inspectors.hpp"

// Test include(s)
#include "detray/test/utils/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GoogleTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <map>

namespace detray {

constexpr std::size_t cache_size{navigation::default_cache_size};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    using context_t = typename navigation_t::detector_type::geometry_context;
    stepping_t _stepping;
    navigation_t _navigation;
    context_t _context{};
};

/// This tests the construction and general methods of the navigator
GTEST_TEST(detray_navigation, direct_navigator_toy_geometry) {
    /*
    using namespace detray;
    using namespace detray::navigation;

    using test_algebra = test::algebra;
    using scalar = test::scalar;
    using point3 = test::point3;
    using vector3 = test::vector3;

    vecmem::host_memory_resource host_mr;

    /// Tolerance for tests
    constexpr double tol{0.01};

    auto [toy_det, names] = build_toy_detector<test_algebra>(host_mr);

    using detector_t = decltype(toy_det);
    using inspector_t = navigation::print_inspector;
    using navigator_t = direct_navigator<detector_t, inspector_t>;
    using constraint_t = constrained_step<scalar>;
    using stepper_t = line_stepper<test_algebra, constraint_t>;

    // Test track
    point3 pos{0.f, 0.f, 0.f};
    vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<test_algebra> traj(pos, 0.f, mom, -1.f);

    navigator_t::sequence_t seq;
    seq.push_back();

    prop_state<stepper_t::state, navigator_t::state> propagation{
        stepper_t::state{traj}, navigator_t::state(toy_det)};
    navigator_t::state &navigation = propagation._navigation;
    stepper_t::state &stepping = propagation._stepping;
    const auto &ctx = propagation._context;
    */

    /*
    stepper_t stepper;
    navigator_t nav;
    navigation::config nav_cfg{};
    stepping::config step_cfg{};

    prop_state<stepper_t::state, navigator_t::state> propagation{
    stepper_t::state{traj}, navigator_t::state(toy_det)};
    navigator_t::state &navigation = propagation._navigation;
    stepper_t::state &stepping = propagation._stepping;
    const auto &ctx = propagation._context;
    */
}

}  // namespace detray