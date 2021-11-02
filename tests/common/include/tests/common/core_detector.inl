/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "detray/core/transform_store.hpp"

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
    detector::mask_container masks;
    geometry::surface_filling_container surfaces = {};

    /// Surface 0
    point3 t0{0., 0., 0.};
    trfs[geometry::e_rectangle2].emplace_back(ctx0, t0);
    masks.template add_mask<geometry::e_rectangle2>(-3., 3.);

    /// Surface 1
    point3 t1{1., 0., 0.};
    trfs[geometry::e_annulus2].emplace_back(ctx0, t1);
    masks.template add_mask<geometry::e_annulus2>(1., 2., 3., 4., 5., 6., 7.);

    /// Surface 2
    point3 t2{2., 0., 0.};
    trfs[geometry::e_trapezoid2].emplace_back(ctx0, t2);
    masks.template add_mask<geometry::e_trapezoid2>(1., 2., 3.);

    detector d("test_detector", host_mr);
    auto &v = d.new_volume({0., 10., -5., 5., -M_PI, M_PI});
    d.add_objects(ctx0, v, surfaces, masks, trfs);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
