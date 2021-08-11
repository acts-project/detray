/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, detector)
{
    using namespace detray;
    using namespace __plugin;

    using detector = detector<>;


    detector d("test_detector");

    //static_transform_store<>::storage static_storage;
    //static_storage.reserve(1);
    static_transform_store<>::context ctx0;

    detector::transform_container static_storages;
    detector::mask_container masks;
    detector::link_container source_links{};
    detector::surface_container surfaces;

    auto &v = d.new_volume({0., 10., -5., 5., -M_PI, M_PI});

    /// Surface 0
    point3 t0{0., 0., 0.};
    std::get<detector::e_rectangle2>(static_storages).emplace_back(t0);
    detector::surface_rectangle rect = {-3., 3.};
    surfaces[detector::e_rectangle2] = {1, detector::e_rectangle2, {0, 1}, 0, dindex_invalid};
    std::get<detector::e_rectangle2>(masks).push_back(rect);
    d.add_surfaces(v, surfaces, masks, static_storages, source_links, ctx0);

    /// Surface 1
    point3 t1{1., 0., 0.};
    std::get<detector::e_annulus2>(static_storages).emplace_back(t1);
    detector::surface_annulus anns = {1., 2., 3., 4., 5., 6., 7.};
    surfaces[detector::e_annulus2] = {1, detector::e_annulus2, {0, 1}, 1, dindex_invalid};
    std::get<detector::e_annulus2>(masks).push_back(anns);
    d.add_surfaces(v, surfaces, masks, static_storages, source_links, ctx0);

    /// Surface 2
    point3 t2{2., 0., 0.};
    std::get<detector::e_trapezoid2>(static_storages).emplace_back(t2);
    detector::surface_trapezoid trap = {1., 2., 3.};
    surfaces[detector::e_trapezoid2] = {1, detector::e_trapezoid2, {0, 1}, 2, dindex_invalid};
    std::get<detector::e_trapezoid2>(masks).push_back(trap);
    d.add_surfaces(v, surfaces, masks, static_storages, source_links, ctx0);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
