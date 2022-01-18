/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/detail/accessor.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    vecmem::host_memory_resource host_mr;
    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;

    auto toy_det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    using context_t = typename decltype(toy_det)::context;
    context_t ctx{};
    auto& volumes = toy_det.volumes();
    auto& surfaces = toy_det.surfaces();
    auto& surfaces_finder = toy_det.get_surfaces_finder();
    auto& transforms = toy_det.transforms();
    auto& masks = toy_det.masks();
    auto& rectangles = masks.template group<0>();
    auto& trapezoids = masks.template group<1>();
    auto& annuli = masks.template group<2>();
    auto& cylinders = masks.template group<3>();
    auto& discs = masks.template group<4>();

    /** source link */
    const dindex inv_sf_finder = dindex_invalid;

    /** Link to outer world (leaving detector) */
    const dindex leaving_world = dindex_invalid;

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 20);
    EXPECT_EQ(surfaces.size(), 3244);
    EXPECT_EQ(surfaces_finder.size(), decltype(toy_det)::N_GRIDS);
    EXPECT_EQ(surfaces_finder.effective_size(),
              n_brl_layers + 2 * n_edc_layers);
    EXPECT_EQ(transforms.size(ctx), 3244);
    EXPECT_EQ(rectangles.size(), 2492);
    EXPECT_EQ(trapezoids.size(), 648);
    EXPECT_EQ(annuli.size(), 0);
    EXPECT_EQ(cylinders.size(), 52);
    EXPECT_EQ(discs.size(), 52);

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
    darray<scalar, 6> bounds = {0., 27., -825., 825., -M_PI, M_PI};
    darray<dindex, 2> range = {0, 16};
    EXPECT_EQ(vol_itr->index(), 0);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of beampipe itself
    range = {0, 1};
    test_module_links(vol_itr->index(), surfaces.begin(), range, range[0],
                      {3, 0}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {1, 8};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 1},
                      {{1, inv_sf_finder},
                       {2, inv_sf_finder},
                       {3, inv_sf_finder},
                       {4, inv_sf_finder},
                       {5, inv_sf_finder},
                       {6, inv_sf_finder},
                       {7, inv_sf_finder}});
    range = {8, 14};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 8},
                      {{14, inv_sf_finder},
                       {15, inv_sf_finder},
                       {16, inv_sf_finder},
                       {17, inv_sf_finder},
                       {18, inv_sf_finder},
                       {19, inv_sf_finder}});

    // disc portals
    range = {14, 16};
    test_portal_links(
        vol_itr->index(), surfaces.begin() + range[0], range, range[0], {4, 0},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -825., -815., -M_PI, M_PI};
    range = {16, 128};
    EXPECT_EQ(vol_itr->index(), 1);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {16, 124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 0}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {124, 126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 14},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {126, 128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 2},
                      {{leaving_world, inv_sf_finder}, {2, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -815., -705., -M_PI, M_PI};
    range = {128, 132};
    EXPECT_EQ(vol_itr->index(), 2);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {128, 130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 16},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {130, 132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 4},
                      {{1, inv_sf_finder}, {3, inv_sf_finder}});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -705., -695., -M_PI, M_PI};
    range = {132, 244};
    EXPECT_EQ(vol_itr->index(), 3);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {132, 240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 108}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {240, 242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 18},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {242, 244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 6},
                      {{2, inv_sf_finder}, {4, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -695., -605., -M_PI, M_PI};
    range = {244, 248};
    EXPECT_EQ(vol_itr->index(), 4);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {244, 246};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 20},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {246, 248};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 8},
                      {{3, inv_sf_finder}, {5, inv_sf_finder}});

    //
    // neg endcap (layer 1)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -605., -595., -M_PI, M_PI};
    range = {248, 360};
    EXPECT_EQ(vol_itr->index(), 5);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {248, 356};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 216}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {356, 358};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 22},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {358, 360};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 10},
                      {{4, inv_sf_finder}, {6, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -595., -500., -M_PI, M_PI};
    range = {360, 370};
    EXPECT_EQ(vol_itr->index(), 6);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {360, 362};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 24},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {362, 370};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 12},
                      {{5, inv_sf_finder},
                       {7, inv_sf_finder},
                       {8, inv_sf_finder},
                       {9, inv_sf_finder},
                       {10, inv_sf_finder},
                       {11, inv_sf_finder},
                       {12, inv_sf_finder},
                       {13, inv_sf_finder}});

    //
    // barrel
    //

    //
    // first layer
    //

    // Check volume
    vol_itr++;
    bounds = {27., 38., -500., 500, -M_PI, M_PI};
    range = {370, 598};
    EXPECT_EQ(vol_itr->index(), 7);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of modules
    range = {370, 594};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 0}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {594, 596};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 26},
                      {{0, inv_sf_finder}, {8, inv_sf_finder}});

    // disc portals
    range = {596, 598};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 20},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {38., 64., -500., 500, -M_PI, M_PI};
    range = {598, 602};
    EXPECT_EQ(vol_itr->index(), 8);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {598, 600};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 28},
                      {{7, inv_sf_finder}, {9, inv_sf_finder}});
    // disc portals
    range = {600, 602};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 22},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // second layer
    //

    // Check volume
    vol_itr++;
    bounds = {64., 80., -500., 500, -M_PI, M_PI};
    range = {602, 1054};
    EXPECT_EQ(vol_itr->index(), 9);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of modules
    range = {602, 1050};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 224}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {1050, 1052};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 30},
                      {{8, inv_sf_finder}, {10, inv_sf_finder}});

    // disc portals
    range = {1052, 1054};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 24},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {80., 108., -500., 500, -M_PI, M_PI};
    range = {1054, 1058};
    EXPECT_EQ(vol_itr->index(), 10);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {1054, 1056};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 32},
                      {{9, inv_sf_finder}, {11, inv_sf_finder}});
    // disc portals
    range = {1056, 1058};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 26},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // third layer
    //

    // Check volume
    vol_itr++;
    bounds = {108., 124., -500., 500, -M_PI, M_PI};
    range = {1058, 1790};
    EXPECT_EQ(vol_itr->index(), 11);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of modules
    range = {1058, 1786};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 672}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {1786, 1788};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 34},
                      {{10, inv_sf_finder}, {12, inv_sf_finder}});

    // disc portals
    range = {1788, 1790};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 28},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {124., 164., -500., 500, -M_PI, M_PI};
    range = {1790, 1794};
    EXPECT_EQ(vol_itr->index(), 12);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {1790, 1792};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 36},
                      {{11, inv_sf_finder}, {13, inv_sf_finder}});
    // disc portals
    range = {1792, 1794};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 30},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // fourth layer
    //

    // Check volume
    vol_itr++;
    bounds = {164., 180., -500., 500, -M_PI, M_PI};
    range = {1794, 2890};
    EXPECT_EQ(vol_itr->index(), 13);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of modules
    range = {1794, 2886};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {0, 1400}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {2886, 2888};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 38},
                      {{12, inv_sf_finder}, {leaving_world, inv_sf_finder}});

    // disc portals
    range = {2888, 2890};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 32},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // positive endcap
    //

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 500., 595., -M_PI, M_PI};
    range = {2890, 2900};
    EXPECT_EQ(vol_itr->index(), 14);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {2890, 2892};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 40},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {2892, 2900};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 34},
                      {{15, inv_sf_finder},
                       {7, inv_sf_finder},
                       {8, inv_sf_finder},
                       {9, inv_sf_finder},
                       {10, inv_sf_finder},
                       {11, inv_sf_finder},
                       {12, inv_sf_finder},
                       {13, inv_sf_finder}});

    //
    // pos endcap (layer 1)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 595., 605., -M_PI, M_PI};
    range = {2900, 3012};
    EXPECT_EQ(vol_itr->index(), 15);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {2900, 3008};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 324}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {3008, 3010};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 42},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3010, 3012};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 42},
                      {{14, inv_sf_finder}, {16, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 605., 695., -M_PI, M_PI};
    range = {3012, 3016};
    EXPECT_EQ(vol_itr->index(), 16);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {3012, 3014};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 44},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3014, 3016};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 44},
                      {{15, inv_sf_finder}, {17, inv_sf_finder}});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 695., 705., -M_PI, M_PI};
    range = {3016, 3128};
    EXPECT_EQ(vol_itr->index(), 17);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {3016, 3124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 432}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {3124, 3126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 46},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3126, 3128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 46},
                      {{16, inv_sf_finder}, {18, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 705., 815., -M_PI, M_PI};
    range = {3128, 3132};
    EXPECT_EQ(vol_itr->index(), 18);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check links of portals
    // cylinder portals
    range = {3128, 3130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 48},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3130, 3132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 48},
                      {{17, inv_sf_finder}, {19, inv_sf_finder}});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 815., 825., -M_PI, M_PI};
    range = {3132, 3244};
    EXPECT_EQ(vol_itr->index(), 19);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);

    // Check the trapezoid modules
    range = {3132, 3240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {1, 540}, {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {3240, 3242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {3, 50},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3242, 3244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {4, 50},
                      {{18, inv_sf_finder}, {leaving_world, inv_sf_finder}});
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
