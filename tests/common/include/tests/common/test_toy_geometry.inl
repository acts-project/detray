/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>

#include "tests/common/read_geometry.hpp"

using namespace detray;

constexpr bool for_surface = true;
constexpr bool for_portal = false;

auto [volumes, surfaces, transforms, cylinders, rectangles] = toy_geometry();

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 4);
    EXPECT_EQ(surfaces.size(), 679);
    EXPECT_EQ(transforms.size(), 679);
    EXPECT_EQ(cylinders.size(), 7);
    EXPECT_EQ(rectangles.size(), 672);

    auto test_portal_links =
        [](dindex vol_index, decltype(surfaces.begin())&& sf_itr,
           darray<dindex, 2>& range, dindex trf_index,
           darray<dindex, 2>&& mask_index, dvector<darray<dindex, 2>>&& edges) {
            for (dindex pti = range[0]; pti < range[1]; pti++) {
                EXPECT_EQ(sf_itr->volume(), vol_index);
                EXPECT_EQ(sf_itr->transform(), trf_index);
                EXPECT_EQ(sf_itr->mask(), mask_index);
                EXPECT_EQ(sf_itr->edge(), edges[pti - range[0]]);
                sf_itr++;
                trf_index++;
                mask_index[1]++;
            }
        };

    auto test_module_links =
        [](dindex vol_index, decltype(surfaces.begin())&& sf_itr,
           darray<dindex, 2>& range, dindex trf_index,
           darray<dindex, 2>&& mask_index, dvector<darray<dindex, 2>>&& edges) {
            for (dindex pti = range[0]; pti < range[1]; pti++) {
                EXPECT_EQ(sf_itr->volume(), vol_index);
                EXPECT_EQ(sf_itr->transform(), trf_index);
                EXPECT_EQ(sf_itr->mask(), mask_index);
                EXPECT_EQ(sf_itr->edge(), edges[0]);
                sf_itr++;
                trf_index++;
                mask_index[1]++;
            }
        };

    //
    // beampipe
    //

    // Check volume
    auto vol_itr = volumes.begin();
    darray<scalar, 6> bounds = {0., 27., -500., 500, -M_PI, M_PI};
    darray<dindex, 2> range = {0, 1};

    EXPECT_EQ(vol_itr->index(), 0);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);

    // Check links of portals
    test_portal_links(vol_itr->index(), surfaces.begin(), range, 0, {0, 0},
                      {{1, dindex_invalid}});

    //
    // first layer
    //

    // Check volume
    vol_itr++;
    bounds = {27., 38., -500., 500, -M_PI, M_PI};
    range = {1, 3};
    EXPECT_EQ(vol_itr->index(), 1);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);

    // Check links of portals
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 1},
                      {{0, dindex_invalid}, {2, dindex_invalid}});

    // Check links of modules
    range = {3, 227};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 0}, {{vol_itr->index(), dindex_invalid}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {38., 64., -500., 500, -M_PI, M_PI};
    range = {227, 229};
    EXPECT_EQ(vol_itr->index(), 2);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);

    // Check links of portals
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 3},
                      {{1, dindex_invalid}, {3, dindex_invalid}});

    //
    // second layer
    //

    // Check volume
    vol_itr++;
    bounds = {64., 80., -500., 500, -M_PI, M_PI};
    range = {229, 231};
    EXPECT_EQ(vol_itr->index(), 3);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);

    // Check links of portals
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 5},
                      {{2, dindex_invalid}, {dindex_invalid, dindex_invalid}});

    // Check links of modules
    range = {231, surfaces.size()};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 224}, {{vol_itr->index(), dindex_invalid}});
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
