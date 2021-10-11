/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "geometry/index_geometry.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, index_geometry) {
    using namespace detray;
    using namespace __plugin;

    using geometry = index_geometry<>;
    using surface = typename geometry::surface;
    using portal = typename geometry::portal;

    geometry g = geometry();

    ASSERT_TRUE(g.n_volumes() == 0);
    ASSERT_TRUE(g.template n_objects<geometry::e_surface>() == 0);
    ASSERT_TRUE(g.template n_objects<geometry::e_portal>() == 0);

    // Add two volumes
    darray<scalar, 6> bounds_0 = {0., 10., -5., 5., -M_PI, M_PI};
    darray<scalar, 6> bounds_1 = {0., 5., -10., 10., -M_PI, M_PI};
    auto &v0 = g.new_volume(bounds_0);
    auto &v1 = g.new_volume(bounds_1);

    ASSERT_TRUE(g.n_volumes() == 2);

    auto &v2 = g.volume_by_index(0);
    auto &v3 = g.volume_by_index(1);

    ASSERT_TRUE(v2.index() == 0);
    ASSERT_TRUE(v3.index() == 1);
    ASSERT_TRUE(v2.bounds() == bounds_0);
    ASSERT_TRUE(v3.bounds() == bounds_1);

    /// volume 0
    /// volume portals
    auto pt0 = portal(0, {geometry::e_portal_cylinder3, {0, 0}}, 0, 1);
    auto pt1 = portal(1, {geometry::e_portal_ring2, {0, 0}}, 0, 2);
    auto pt2 = portal(2, {geometry::e_portal_ring2, {1, 1}}, 0, 3);
    /// Surface 0
    auto sf0 = surface(3, {geometry::e_rectangle2, 0}, 0, 4);
    /// Surface 1
    auto sf1 = surface(4, {geometry::e_trapezoid2, 0}, 0, 5);

    /// volume 1
    /// volume portals
    auto pt3 = portal(0, {geometry::e_portal_cylinder3, {0, 1}}, 1, 6);
    auto pt4 = portal(1, {geometry::e_portal_ring2, {0, 1}}, 1, 7);
    /// Surface 2
    auto sf2 = surface(0, {geometry::e_rectangle2, 0}, 1, 8);
    /// Surface 3
    auto sf3 = surface(1, {geometry::e_trapezoid2, 0}, 1, 9);

    /// Correct links for volume 1
    dindex trf_offset_vol1 = 5;
    dindex mask_offset_cyl = 1;
    dindex mask_offset_rg = 2;
    dindex mask_offset_rect = 1;
    dindex mask_offset_trap = 1;

    // update transform links
    g.update_transform_link(pt3, trf_offset_vol1);
    g.update_transform_link(pt4, trf_offset_vol1++);
    g.update_transform_link(sf2, ++trf_offset_vol1);
    g.update_transform_link(sf3, trf_offset_vol1);

    ASSERT_EQ(pt3.transform(), 5);
    ASSERT_EQ(pt4.transform(), 6);
    ASSERT_EQ(sf2.transform(), 7);
    ASSERT_EQ(sf3.transform(), 8);

    // update mask links
    g.update_mask_link(pt3, mask_offset_cyl);
    g.update_mask_link(pt4, mask_offset_rg);
    g.update_mask_link(sf2, mask_offset_rect);
    g.update_mask_link(sf3, mask_offset_trap);

    ASSERT_EQ(std::get<1>(pt3.mask())[0], 1);
    ASSERT_EQ(std::get<1>(pt3.mask())[1], 2);
    ASSERT_EQ(std::get<1>(pt4.mask())[0], 2);
    ASSERT_EQ(std::get<1>(pt4.mask())[1], 3);
    ASSERT_EQ(std::get<1>(sf2.mask()), 1);
    ASSERT_EQ(std::get<1>(sf3.mask()), 1);

    /// fill surfaces and portals (through volume alias names)
    dvector<portal> portals_vol0 = {pt0, pt1, pt2};
    dvector<portal> portals_vol1 = {pt3, pt4};
    dvector<surface> surfaces_vol0 = {sf0, sf1};
    dvector<surface> surfaces_vol1 = {sf2, sf3};

    g.template add_objects<geometry::e_portal>(v2, portals_vol0);
    g.template add_objects<geometry::e_surface>(v2, surfaces_vol0);
    g.template add_objects<geometry::e_portal>(v3, portals_vol1);
    g.template add_objects<geometry::e_surface>(v3, surfaces_vol1);

    // Are the surfaces/portals filled correctly?
    ASSERT_EQ(g.template n_objects<geometry::e_portal>(), 5);
    ASSERT_EQ(g.template n_objects<geometry::e_surface>(), 4);

    // Are the ranges updated correctly?
    auto objects_range = darray<dindex, 2>{0, 3};
    ASSERT_TRUE(v2.template range<geometry::e_portal>() == objects_range);
    objects_range = darray<dindex, 2>{0, 2};
    ASSERT_TRUE(v2.template range<geometry::e_surface>() == objects_range);
    objects_range = darray<dindex, 2>{3, 5};
    ASSERT_TRUE(v3.template range<geometry::e_portal>() == objects_range);
    objects_range = darray<dindex, 2>{2, 4};
    ASSERT_TRUE(v3.template range<geometry::e_surface>() == objects_range);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
