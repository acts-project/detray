/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray include(s)
#include "detray/propagator/line_stepper.hpp"

#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s)
#include <memory>

// google-test include(s)
#include <gtest/gtest.h>

using namespace detray;

using vector3 = test::vector3;
using point3 = test::point3;
using transform3 = test::transform3;

namespace {

constexpr scalar tol{1e-3f};

vecmem::host_memory_resource host_mr;

// dummy navigation struct
struct nav_state {
    /// New detector
    nav_state(vecmem::host_memory_resource &mr)
        : m_step_size{1.f * unit<scalar>::mm},
          m_det{std::make_unique<detray::detector<>>(mr)} {
        // Empty dummy volume
        volume_builder<detray::detector<>> vbuilder{volume_id::e_cylinder};
        vbuilder.build(*m_det);
    }

    scalar operator()() const { return m_step_size; }
    inline auto current_object() const -> dindex { return dindex_invalid; }
    inline auto tolerance() const -> scalar { return tol; }
    inline auto detector() const -> const detray::detector<> * {
        return m_det.get();
    }
    inline auto volume() -> unsigned int { return 0u; }
    inline void set_full_trust() {}
    inline void set_high_trust() {}
    inline void set_fair_trust() {}
    inline void set_no_trust() {}
    inline bool abort() { return false; }

    scalar m_step_size;
    std::unique_ptr<detray::detector<>> m_det;
};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    stepping_t _stepping;
    navigation_t _navigation;
};

}  // namespace

// This tests the base functionality of the line stepper
GTEST_TEST(detray_propagator, line_stepper) {
    using namespace step;

    // Line stepper with and without constrained stepping
    using line_stepper_t = line_stepper<transform3>;
    using cline_stepper_t = line_stepper<transform3, constrained_step<>>;

    point3 pos{0.f, 0.f, 0.f};
    vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<transform3> track(pos, 0.f, mom, -1.f);
    free_track_parameters<transform3> c_track(pos, 0.f, mom, -1.f);

    line_stepper_t l_stepper;
    cline_stepper_t cl_stepper;

    prop_state<line_stepper_t::state, nav_state> propagation{
        line_stepper_t::state{track}, nav_state{host_mr}};
    prop_state<cline_stepper_t::state, nav_state> c_propagation{
        cline_stepper_t::state{c_track}, nav_state{host_mr}};

    cline_stepper_t::state &cl_state = c_propagation._stepping;

    // Test the setting of step constraints
    cl_state.template set_constraint<constraint::e_accuracy>(10.f *
                                                             unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_actor>(2.f *
                                                          unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_aborter>(5.f *
                                                            unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_user>(0.5f *
                                                         unit<scalar>::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5f * unit<scalar>::mm, tol);

    // Release all except 'actor', then set 'user' again
    cl_state.template release_step<constraint::e_accuracy>();
    cl_state.template release_step<constraint::e_aborter>();
    cl_state.template release_step<constraint::e_user>();
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                2.f * unit<scalar>::mm, tol);

    cl_state.template set_constraint<constraint::e_user>(0.5f *
                                                         unit<scalar>::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5f * unit<scalar>::mm, tol);

    // Run a few steps
    ASSERT_TRUE(l_stepper.step(propagation));
    // Step constraint to half step size
    ASSERT_TRUE(cl_stepper.step(c_propagation));
    ASSERT_TRUE(cl_stepper.step(c_propagation));

    track = propagation._stepping();
    ASSERT_NEAR(track.pos()[0], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(track.pos()[1], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(track.pos()[2], 0.f, tol);

    c_track = c_propagation._stepping();
    ASSERT_NEAR(c_track.pos()[0], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(c_track.pos()[1], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(c_track.pos()[2], 0.f, tol);

    ASSERT_TRUE(l_stepper.step(propagation));

    track = propagation._stepping();
    ASSERT_NEAR(track.pos()[0], constant<scalar>::sqrt2, tol);
    ASSERT_NEAR(track.pos()[1], constant<scalar>::sqrt2, tol);
    ASSERT_NEAR(track.pos()[2], 0.f, tol);
}
