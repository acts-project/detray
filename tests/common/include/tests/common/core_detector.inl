/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "core/detector.hpp"
#include "core/transform_store.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, detector) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    using detector = detector<>;
    using geometry = detector::geometry;

    static_transform_store<>::context ctx0;

    detector::transform_container trfs;
    geometry::surface_mask_container masks;
    geometry::surface_filling_container surfaces = {};

    /// Surface 0
    point3 t0{0., 0., 0.};
    trfs[geometry::surface_rectangle::mask_context].emplace_back(ctx0, t0);
    geometry::surface_rectangle rect = {-3., 3.};
    std::get<geometry::surface_rectangle::mask_context>(masks).push_back(rect);

    /// Surface 1
    point3 t1{1., 0., 0.};
    trfs[geometry::surface_annulus::mask_context].emplace_back(ctx0, t1);
    geometry::surface_annulus anns = {1., 2., 3., 4., 5., 6., 7.};
    std::get<geometry::surface_annulus::mask_context>(masks).push_back(anns);

    /// Surface 2
    point3 t2{2., 0., 0.};
    trfs[geometry::surface_trapezoid::mask_context].emplace_back(ctx0, t2);
    geometry::surface_trapezoid trap = {1., 2., 3.};
    std::get<geometry::surface_trapezoid::mask_context>(masks).push_back(trap);

    detector d("test_detector", host_mr);
    auto &v = d.new_volume("test_volume", {0., 10., -5., 5., -M_PI, M_PI});
    d.template add_objects<geometry::e_surface>(v, surfaces, masks, trfs, ctx0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
