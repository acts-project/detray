/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "tests/common/read_geometry.hpp"

using namespace detray;

auto [volumes, surfaces, transforms, cylinders, discs, rectangles] =
    toy_geometry();

constexpr auto for_surface = decltype(volumes)::value_type::objects::e_surface;
constexpr auto for_portal = decltype(volumes)::value_type::objects::e_portal;

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    /** source link */
    const dindex inv_sf_finder = dindex_invalid;

    /** Link to outer world (leaving detector) */
    const dindex leaving_world = dindex_invalid;

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 4);
    EXPECT_EQ(surfaces.size(), 687);
    EXPECT_EQ(transforms.size(), 687);
    EXPECT_EQ(cylinders.size(), 7);
    EXPECT_EQ(discs.size(), 8);
    EXPECT_EQ(rectangles.size(), 672);

    /** Test the links of portals (into the next volume or invalid if we leave
     * the detector).
     *
     * @param vol_index volume the portals belong to
     * @param sf_itr iterator into the surface container, start of the portals
     * @param range index range of the portals in the surface container
     * @param trf_index index of the transform (trf container) for the portal
     * @param mask_index type and index of portal mask in respective mask cont
     * @param edges links to next volume and next surfaces finder
     */
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

    /** Test the links of module surface (alway stay in their volume).
     *
     * @param vol_index volume the modules belong to
     * @param sf_itr iterator into the surface container, start of the modules
     * @param range index range of the modules in the surface container
     * @param trf_index index of the transform (trf container) for the module
     * @param mask_index type and index of module mask in respective mask cont
     * @param edges links to next volume and next surfaces finder
     */
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
    darray<dindex, 2> range = {0, 3};
    EXPECT_EQ(vol_itr->index(), 0);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);
    range = {dindex_invalid, dindex_invalid};
    EXPECT_EQ(vol_itr->template range<for_surface>(), range);

    // Check links of portals
    // disc portals
    range = {0, 2};
    test_portal_links(
        vol_itr->index(), surfaces.begin(), range, range[0], {0, 0},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // cylinder portals
    range = {2, 3};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 0}, {{1, inv_sf_finder}});

    //
    // first layer
    //

    // Check volume
    vol_itr++;
    bounds = {27., 38., -500., 500, -M_PI, M_PI};
    range = {3, 7};
    EXPECT_EQ(vol_itr->index(), 1);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);
    range = {7, 231};
    EXPECT_EQ(vol_itr->template range<for_surface>(), range);

    // Check links of portals
    // disc portals
    range = {3, 5};
    test_portal_links(
        vol_itr->index(), surfaces.begin() + range[0], range, range[0], {0, 2},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // cylinder portals
    range = {5, 7};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 1},
                      {{0, inv_sf_finder}, {2, inv_sf_finder}});

    // Check links of modules
    range = {7, 231};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {2, 0}, {{vol_itr->index(), inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {38., 64., -500., 500, -M_PI, M_PI};
    range = {231, 235};
    EXPECT_EQ(vol_itr->index(), 2);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);
    range = {dindex_invalid, dindex_invalid};
    EXPECT_EQ(vol_itr->template range<for_surface>(), range);

    // Check links of portals
    // disc portals
    range = {231, 233};
    test_portal_links(
        vol_itr->index(), surfaces.begin() + range[0], range, range[0], {0, 4},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // cylinder portals
    range = {233, 235};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 3},
                      {{1, inv_sf_finder}, {3, inv_sf_finder}});

    //
    // second layer
    //

    // Check volume
    vol_itr++;
    bounds = {64., 80., -500., 500, -M_PI, M_PI};
    range = {235, 239};
    EXPECT_EQ(vol_itr->index(), 3);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->template range<for_portal>(), range);
    range = {239, surfaces.size()};
    EXPECT_EQ(vol_itr->template range<for_surface>(), range);

    // Check links of portals
    // disc portals
    range = {235, 237};
    test_portal_links(
        vol_itr->index(), surfaces.begin() + range[0], range, range[0], {0, 6},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // cylinder portals
    range = {237, 239};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 5},
                      {{2, inv_sf_finder}, {leaving_world, inv_sf_finder}});

    // Check links of modules
    range = {239, surfaces.size()};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {2, 224}, {{vol_itr->index(), inv_sf_finder}});
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
