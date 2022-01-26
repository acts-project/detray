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

    using detector_t = decltype(toy_det);
    using volume_t = typename detector_t::volume_type;
    using mask_defs = typename detector_t::mask_defs;
    using mask_link_t = typename mask_defs::link_type;
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
    /** Mask ids */
    constexpr auto cylinder_id = mask_defs::e_portal_cylinder3;
    constexpr auto disc_id = mask_defs::e_portal_ring2;
    constexpr auto rectangle_id = mask_defs::e_rectangle2;
    constexpr auto trapezoid_id = mask_defs::e_trapezoid2;

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
           darray<dindex, 2>& range, dindex trf_index, mask_link_t&& mask_index,
           dvector<darray<dindex, 2>>&& edges) {
            for (dindex pti = range[0]; pti < range[1]; ++pti) {
                EXPECT_EQ(sf_itr->volume(), vol_index);
                EXPECT_EQ(sf_itr->transform(), trf_index);
                EXPECT_EQ(sf_itr->mask(), mask_index);
                EXPECT_EQ(sf_itr->edge(), edges[pti - range[0]]);
                ++sf_itr;
                ++trf_index;
                ++std::get<1>(mask_index);
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
           darray<dindex, 2>& range, dindex trf_index, mask_link_t&& mask_index,
           dvector<darray<dindex, 2>>&& edges) {
            for (dindex pti = range[0]; pti < range[1]; ++pti) {
                EXPECT_EQ(sf_itr->volume(), vol_index);
                EXPECT_EQ(sf_itr->transform(), trf_index);
                EXPECT_EQ(sf_itr->mask(), mask_index);
                EXPECT_EQ(sf_itr->edge(), edges[0]);
                ++sf_itr;
                ++trf_index;
                ++std::get<1>(mask_index);
            }
        };

    /** Test the surface grid.
     *
     * @param volume volume of detector being tested
     * @param sf_finder surface finder of detector
     * @param sf_container surface container of detector
     * @param range index range of the modules in the surface container
     */
    auto test_surfaces_grid =
        [](decltype(volumes.begin())& vol_itr,
           const typename detector_t::surfaces_finder_type& sf_finder,
           const typename detector_t::surface_container& sf_container,
           darray<dindex, 2>& range) {
            auto sf_grid = sf_finder[vol_itr->sf_finder_index()];

            std::vector<dindex> indices;

            // Check if right volume is linked to surface in grid
            for (unsigned int i = 0; i < sf_grid.axis_p0().bins(); i++) {
                for (unsigned int j = 0; j < sf_grid.axis_p1().bins(); j++) {
                    auto sf_indices = sf_grid.bin(i, j);
                    for (auto sf_idx : sf_indices) {
                        EXPECT_EQ(sf_container[sf_idx].volume(),
                                  vol_itr->index());
                        indices.push_back(sf_idx);
                    }
                }
            }

            // Check if surface indices are consistent with given range of
            // modules
            EXPECT_EQ(indices.size(), range[1] - range[0]);
            for (dindex pti = range[0]; pti < range[1]; pti++) {
                EXPECT_TRUE(std::find(indices.begin(), indices.end(), pti) !=
                            indices.end());
            }
        };

    //
    // beampipe
    //

    // Check volume
    auto vol_itr = volumes.begin();
    darray<scalar, 6> bounds = {0., 27., -825., 825., -M_PI, M_PI};
    darray<dindex, 2> range = {0, 16};
    volume_t::sf_finder_link_t sf_finder_link = {
        volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 0);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of beampipe itself
    range = {0, 1};
    test_module_links(vol_itr->index(), surfaces.begin(), range, range[0],
                      {mask_defs::e_cylinder3, 0},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check links of portals
    // cylinder portals
    range = {1, 8};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 1},
                      {{1, inv_sf_finder},
                       {2, inv_sf_finder},
                       {3, inv_sf_finder},
                       {4, inv_sf_finder},
                       {5, inv_sf_finder},
                       {6, inv_sf_finder},
                       {7, inv_sf_finder}});
    range = {8, 14};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 8},
                      {{14, inv_sf_finder},
                       {15, inv_sf_finder},
                       {16, inv_sf_finder},
                       {17, inv_sf_finder},
                       {18, inv_sf_finder},
                       {19, inv_sf_finder}});

    // disc portals
    range = {14, 16};
    test_portal_links(
        vol_itr->index(), surfaces.begin() + range[0], range, range[0],
        {disc_id, 0},
        {{leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -825., -815., -M_PI, M_PI};
    range = {16, 128};
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 0};
    EXPECT_EQ(vol_itr->index(), 1);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {16, 124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 0},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {124, 126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 14},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {126, 128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 2},
                      {{leaving_world, inv_sf_finder}, {2, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -815., -705., -M_PI, M_PI};
    range = {128, 132};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 2);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {128, 130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 16},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {130, 132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 4},
                      {{1, inv_sf_finder}, {3, inv_sf_finder}});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -705., -695., -M_PI, M_PI};
    range = {132, 244};
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 1};
    EXPECT_EQ(vol_itr->index(), 3);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {132, 240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 108},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {240, 242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 18},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {242, 244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 6},
                      {{2, inv_sf_finder}, {4, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -695., -605., -M_PI, M_PI};
    range = {244, 248};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 4);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {244, 246};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 20},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {246, 248};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 8},
                      {{3, inv_sf_finder}, {5, inv_sf_finder}});

    //
    // neg endcap (layer 1)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -605., -595., -M_PI, M_PI};
    range = {248, 360};
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 2};
    EXPECT_EQ(vol_itr->index(), 5);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {248, 356};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 216},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {356, 358};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 22},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {358, 360};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 10},
                      {{4, inv_sf_finder}, {6, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., -595., -500., -M_PI, M_PI};
    range = {360, 370};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 6);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {360, 362};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 24},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {362, 370};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 12},
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
    sf_finder_link = {volume_t::sf_finders::e_z_phi_grid, 3};
    EXPECT_EQ(vol_itr->index(), 7);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of modules
    range = {370, 594};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {rectangle_id, 0},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {594, 596};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 26},
                      {{0, inv_sf_finder}, {8, inv_sf_finder}});

    // disc portals
    range = {596, 598};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 20},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {38., 64., -500., 500, -M_PI, M_PI};
    range = {598, 602};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 8);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);
    // Check links of portals
    // cylinder portals
    range = {598, 600};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 28},
                      {{7, inv_sf_finder}, {9, inv_sf_finder}});
    // disc portals
    range = {600, 602};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 22},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // second layer
    //

    // Check volume
    vol_itr++;
    bounds = {64., 80., -500., 500, -M_PI, M_PI};
    range = {602, 1054};
    sf_finder_link = {volume_t::sf_finders::e_z_phi_grid, 4};
    EXPECT_EQ(vol_itr->index(), 9);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of modules
    range = {602, 1050};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {rectangle_id, 224},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1050, 1052};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 30},
                      {{8, inv_sf_finder}, {10, inv_sf_finder}});

    // disc portals
    range = {1052, 1054};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 24},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {80., 108., -500., 500, -M_PI, M_PI};
    range = {1054, 1058};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 10);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);
    // Check links of portals
    // cylinder portals
    range = {1054, 1056};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 32},
                      {{9, inv_sf_finder}, {11, inv_sf_finder}});
    // disc portals
    range = {1056, 1058};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 26},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // third layer
    //

    // Check volume
    vol_itr++;
    bounds = {108., 124., -500., 500, -M_PI, M_PI};
    range = {1058, 1790};
    sf_finder_link = {volume_t::sf_finders::e_z_phi_grid, 5};
    EXPECT_EQ(vol_itr->index(), 11);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of modules
    range = {1058, 1786};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {rectangle_id, 672},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1786, 1788};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 34},
                      {{10, inv_sf_finder}, {12, inv_sf_finder}});

    // disc portals
    range = {1788, 1790};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 28},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {124., 164., -500., 500, -M_PI, M_PI};
    range = {1790, 1794};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 12);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {1790, 1792};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 36},
                      {{11, inv_sf_finder}, {13, inv_sf_finder}});
    // disc portals
    range = {1792, 1794};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 30},
                      {{6, inv_sf_finder}, {14, inv_sf_finder}});

    //
    // fourth layer
    //

    // Check volume
    vol_itr++;
    bounds = {164., 180., -500., 500, -M_PI, M_PI};
    range = {1794, 2890};
    sf_finder_link = {volume_t::sf_finders::e_z_phi_grid, 6};
    EXPECT_EQ(vol_itr->index(), 13);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of modules
    range = {1794, 2886};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {rectangle_id, 1400},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {2886, 2888};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 38},
                      {{12, inv_sf_finder}, {leaving_world, inv_sf_finder}});

    // disc portals
    range = {2888, 2890};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 32},
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
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 14);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {2890, 2892};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 40},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {2892, 2900};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 34},
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
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 7};
    EXPECT_EQ(vol_itr->index(), 15);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {2900, 3008};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 324},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3008, 3010};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 42},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3010, 3012};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 42},
                      {{14, inv_sf_finder}, {16, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 605., 695., -M_PI, M_PI};
    range = {3012, 3016};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 16);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3012, 3014};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 44},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3014, 3016};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 44},
                      {{15, inv_sf_finder}, {17, inv_sf_finder}});

    //
    // pos endcap (layer 2)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 695., 705., -M_PI, M_PI};
    range = {3016, 3128};
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 8};
    EXPECT_EQ(vol_itr->index(), 17);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {3016, 3124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 432},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3124, 3126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 46},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3126, 3128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 46},
                      {{16, inv_sf_finder}, {18, inv_sf_finder}});

    //
    // gap
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 705., 815., -M_PI, M_PI};
    range = {3128, 3132};
    sf_finder_link = {volume_t::sf_finders::e_unknown, inv_sf_finder};
    EXPECT_EQ(vol_itr->index(), 18);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3128, 3130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 48},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3130, 3132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 48},
                      {{17, inv_sf_finder}, {19, inv_sf_finder}});

    //
    // pos endcap (layer 3)
    //

    // Check volume
    vol_itr++;
    bounds = {27., 180., 815., 825., -M_PI, M_PI};
    range = {3132, 3244};
    sf_finder_link = {volume_t::sf_finders::e_r_phi_grid, 9};
    EXPECT_EQ(vol_itr->index(), 19);
    EXPECT_EQ(vol_itr->bounds(), bounds);
    EXPECT_EQ(vol_itr->range(), range);
    EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);

    // Check the trapezoid modules
    range = {3132, 3240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {trapezoid_id, 540},
                      {{vol_itr->index(), inv_sf_finder}});

    // Check link of surfaces in surface finder
    test_surfaces_grid(vol_itr, surfaces_finder, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3240, 3242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {cylinder_id, 50},
                      {{0, inv_sf_finder}, {leaving_world, inv_sf_finder}});
    // disc portals
    range = {3242, 3244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {disc_id, 50},
                      {{18, inv_sf_finder}, {leaving_world, inv_sf_finder}});
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
