/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#include <gtest/gtest.h>

#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

namespace {

/// Functor that returns the volume/sf_finder links of a particluar mask
/// instance in the mask container of a detector
struct volume_link_getter {
    using output_type = dindex;

    template <typename mask_group_t, typename mask_range_t>
    inline output_type operator()(mask_group_t& mask_group,
                                  mask_range_t& mask_range) {
        return mask_group[mask_range].volume_link();
    }
};

/// Functor that tests the volume indices of surfaces in the grid
struct surface_grid_tester {
    using output_type = bool;

    /// Call for grid types
    template <
        typename grid_group_t, typename surface_container_t,
        std::enable_if_t<not std::is_same_v<typename grid_group_t::value_type,
                                            brute_force_finder>,
                         bool> = true>
    inline output_type operator()(const grid_group_t& grid_group,
                                  const dindex grid_index,
                                  const dindex volume_index,
                                  const surface_container_t& surfaces,
                                  const darray<dindex, 2>& /*range*/) {

        const auto& sf_grid = grid_group.at(grid_index);

        std::vector<dindex> indices;

        // Check if right volume is linked to surface in grid
        for (unsigned int i = 0; i < sf_grid.axis_p0().bins(); ++i) {
            for (unsigned int j = 0; j < sf_grid.axis_p1().bins(); ++j) {
                auto sf_indices = sf_grid.bin(i, j);
                for (const auto sf_idx : sf_indices) {
                    EXPECT_EQ(surfaces.at(sf_idx).volume(), volume_index);
                    indices.push_back(sf_idx);
                }
            }
        }

        // Check if surface indices are consistent with given range of
        // modules
        /*EXPECT_EQ(indices.size(), range[1] - range[0]);
        for (dindex pti = range[0]; pti < range[1]; ++pti) {
            EXPECT_TRUE(std::find(indices.begin(), indices.end(), pti) !=
                        indices.end());
        }*/

        return true;
    }

    /// Call for the brute force type (do nothing)
    template <typename grid_group_t, typename surface_container_t,
              std::enable_if_t<std::is_same_v<typename grid_group_t::value_type,
                                              brute_force_finder>,
                               bool> = true>
    inline output_type operator()(const grid_group_t&, const dindex,
                                  const dindex, const surface_container_t&,
                                  const darray<dindex, 2>&) {
        return true;
    }
};

}  // anonymous namespace

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    vecmem::host_memory_resource host_mr;
    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;

    const auto toy_det =
        create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    using detector_t = decltype(toy_det);
    using volume_t = typename detector_t::volume_type;
    using context_t = typename decltype(toy_det)::context;
    using mask_ids = typename detector_t::masks::id;
    using mask_link_t = typename detector_t::surface_type::mask_link;
    using material_link_t = typename detector_t::surface_type::material_link;
    using material_ids = typename detector_t::materials::id;
    using sf_finder_ids = typename detector_t::sf_finders::id;
    using sf_finder_link_t = typename volume_t::sf_finder_link_type;

    context_t ctx{};
    auto& volumes = toy_det.volumes();
    auto& surfaces = toy_det.surfaces();
    // auto& sf_finders = toy_det.sf_finder_store();
    auto& transforms = toy_det.transform_store();
    auto& masks = toy_det.mask_store();
    auto& materials = toy_det.material_store();

    // Materials
    auto portal_mat =
        material_slab<scalar>(vacuum<scalar>(), 0. * unit_constants::mm);
    auto beampipe_mat = material_slab<scalar>(beryllium_tml<scalar>(),
                                              0.8 * unit_constants::mm);
    auto pixel_mat =
        material_slab<scalar>(silicon_tml<scalar>(), 0.15 * unit_constants::mm);

    /** Link to outer world (leaving detector) */
    const dindex leaving_world = dindex_invalid;

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 20);
    EXPECT_EQ(surfaces.size(), 3244);
    /*EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_brute_force>(), 1);
    EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_z_phi_grid>(),
              n_brl_layers);
    EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_r_phi_grid>(), 14);*/
    EXPECT_EQ(transforms.size(ctx), 3244);
    EXPECT_EQ(masks.template size<mask_ids::e_rectangle2>(), 2492);
    EXPECT_EQ(masks.template size<mask_ids::e_trapezoid2>(), 648);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_cylinder3>(), 52);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_ring2>(), 52);
    EXPECT_EQ(materials.template size<material_ids::e_slab>(), 3244);

    /** Test the links of a volume.
     *
     * @param vol_index volume the modules belong to
     * @param sf_itr iterator into the surface container, start of the modules
     * @param range index range of the modules in the surface container
     * @param trf_index index of the transform (trf container) for the module
     * @param mask_index type and index of module mask in respective mask cont
     * @param volume_links links to next volume and next surfaces finder
     */
    auto test_volume_links =
        [&](decltype(volumes.begin())& vol_itr, const dindex vol_index,
            const darray<scalar, 6>& bounds, const darray<dindex, 2>& range,
            const sf_finder_link_t& /*sf_finder_link*/) {
            EXPECT_EQ(vol_itr->index(), vol_index);
            EXPECT_EQ(vol_itr->bounds(), bounds);
            EXPECT_EQ(vol_itr->range(), range);
            // EXPECT_EQ(vol_itr->sf_finder_link(), sf_finder_link);
        };

    /** Test the links of portals (into the next volume or invalid if we leave
     * the detector).
     *
     * @param vol_index volume the portals belong to
     * @param sf_itr iterator into the surface container, start of the portals
     * @param range index range of the portals in the surface container
     * @param trf_index index of the transform (trf container) for the portal
     * @param mask_index type and index of portal mask in respective mask cont
     * @param volume_links links to next volume contained in the masks
     */
    auto test_portal_links = [&](const dindex vol_index,
                                 decltype(surfaces.begin())&& sf_itr,
                                 const darray<dindex, 2>& range,
                                 dindex trf_index, mask_link_t&& mask_link,
                                 material_link_t&& material_index,
                                 const material_slab<scalar>& mat,
                                 const dvector<dindex>&& volume_links) {
        for (dindex pti = range[0]; pti < range[1]; ++pti) {
            EXPECT_EQ(sf_itr->volume(), vol_index);
            EXPECT_EQ(sf_itr->transform(), trf_index);
            EXPECT_EQ(sf_itr->mask(), mask_link);
            const auto volume_link =
                masks.template call<volume_link_getter>(sf_itr->mask());
            EXPECT_EQ(volume_link, volume_links[pti - range[0]]);
            EXPECT_EQ(
                materials
                    .group<material_ids::e_slab>()[sf_itr->material_range()],
                mat);

            ++sf_itr;
            ++trf_index;
            ++mask_link;
            ++material_index;
        }
    };

    /** Test the links of module surface (alway stay in their volume).
     *
     * @param vol_index volume the modules belong to
     * @param sf_itr iterator into the surface container, start of the modules
     * @param range index range of the modules in the surface container
     * @param trf_index index of the transform (trf container) for the module
     * @param mask_index type and index of module mask in respective mask cont
     * @param volume_links links to next volume contained in the masks
     */
    auto test_module_links = [&](const dindex vol_index,
                                 decltype(surfaces.begin())&& sf_itr,
                                 const darray<dindex, 2>& range,
                                 dindex trf_index, mask_link_t&& mask_index,
                                 material_link_t&& material_index,
                                 const material_slab<scalar>& mat,
                                 const dvector<dindex>&& volume_links) {
        for (dindex pti = range[0]; pti < range[1]; ++pti) {
            EXPECT_EQ(sf_itr->volume(), vol_index);
            EXPECT_EQ(sf_itr->transform(), trf_index);
            EXPECT_EQ(sf_itr->mask(), mask_index);
            EXPECT_EQ(sf_itr->material(), material_index);
            const auto volume_link =
                masks.template call<volume_link_getter>(sf_itr->mask());
            EXPECT_EQ(volume_link, volume_links[0]);
            EXPECT_EQ(
                materials
                    .group<material_ids::e_slab>()[sf_itr->material_range()],
                mat);
            ++sf_itr;
            ++trf_index;
            ++mask_index;
            ++material_index;
        }
    };

    /** Test the surface grid.
     *
     * @param volume volume of detector being tested
     * @param sf_finders surface finder container of detector
     * @param surfaces surface container of detector
     * @param range index range of the modules in the surface container
     */
    /*auto test_surfaces_grid =
        [](decltype(volumes.begin())& vol_itr,
           const typename detector_t::sf_finder_container& sf_finder_cont,
           const typename detector_t::surface_container& surface_cont,
           darray<dindex, 2>& range) {
            // Call test functor
            sf_finder_cont.template call<surface_grid_tester>(
                vol_itr->sf_finder_link(), vol_itr->index(), surface_cont,
                range);
        };*/

    //
    // beampipe
    //

    // Check volume
    auto vol_itr = volumes.begin();
    darray<scalar, 6> bounds = {0., 27., -825., 825., -M_PI, M_PI};
    darray<dindex, 2> range = {0, 16};
    sf_finder_link_t sf_finder_link{sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 0, bounds, range, sf_finder_link);

    // Check links of beampipe itself
    range = {0, 1};
    test_module_links(vol_itr->index(), surfaces.begin(), range, range[0],
                      {mask_ids::e_cylinder3, 0}, {material_ids::e_slab, 0},
                      beampipe_mat, {vol_itr->index()});

    // Check links of portals
    // cylinder portals
    range = {1, 8};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 1},
                      {material_ids::e_slab, 1}, portal_mat,
                      {1, 2, 3, 4, 5, 6, 7});
    range = {8, 14};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 8},
                      {material_ids::e_slab, 8}, portal_mat,
                      {14, 15, 16, 17, 18, 19});

    // disc portals
    range = {14, 16};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 0},
                      {material_ids::e_slab, 14}, portal_mat,
                      {leaving_world, leaving_world});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -825., -815., -M_PI, M_PI};
    range = {16, 128};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 1, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {16, 124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 0},
                      {material_ids::e_slab, 16}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {124, 126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 14},
                      {material_ids::e_slab, 124}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {126, 128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 2},
                      {material_ids::e_slab, 126}, portal_mat,
                      {leaving_world, 2});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -815., -705., -M_PI, M_PI};
    range = {128, 132};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 2, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {128, 130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 16},
                      {material_ids::e_slab, 128}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {130, 132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 4},
                      {material_ids::e_slab, 130}, portal_mat, {1, 3});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -705., -695., -M_PI, M_PI};
    range = {132, 244};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 1};

    // Test the links in the volumes
    test_volume_links(vol_itr, 3, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {132, 240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 108},
                      {material_ids::e_slab, 132}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {240, 242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 18},
                      {material_ids::e_slab, 240}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {242, 244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 6},
                      {material_ids::e_slab, 242}, portal_mat, {2, 4});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -695., -605., -M_PI, M_PI};
    range = {244, 248};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 4, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {244, 246};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 20},
                      {material_ids::e_slab, 244}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {246, 248};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 8},
                      {material_ids::e_slab, 246}, portal_mat, {3, 5});

    //
    // neg endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -605., -595., -M_PI, M_PI};
    range = {248, 360};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 2};

    // Test the links in the volumes
    test_volume_links(vol_itr, 5, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {248, 356};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 216},
                      {material_ids::e_slab, 248}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {356, 358};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 22},
                      {material_ids::e_slab, 356}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {358, 360};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 10},
                      {material_ids::e_slab, 358}, portal_mat, {4, 6});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., -595., -500., -M_PI, M_PI};
    range = {360, 370};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 6, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {360, 362};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 24},
                      {material_ids::e_slab, 360}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {362, 370};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 12},
                      {material_ids::e_slab, 362}, portal_mat,
                      {5, 7, 8, 9, 10, 11, 12, 13});

    //
    // barrel
    //

    //
    // first layer
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 38., -500., 500, -M_PI, M_PI};
    range = {370, 598};
    sf_finder_link = {sf_finder_ids::e_z_phi_grid, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 7, bounds, range, sf_finder_link);

    // Check links of modules
    range = {370, 594};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_rectangle2, 0},
                      {material_ids::e_slab, 370}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {594, 596};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 26},
                      {material_ids::e_slab, 594}, portal_mat, {0, 8});

    // disc portals
    range = {596, 598};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 20},
                      {material_ids::e_slab, 596}, portal_mat, {6, 14});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {38., 64., -500., 500, -M_PI, M_PI};
    range = {598, 602};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 8, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {598, 600};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 28},
                      {material_ids::e_slab, 598}, portal_mat, {7, 9});
    // disc portals
    range = {600, 602};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 22},
                      {material_ids::e_slab, 600}, portal_mat, {6, 14});

    //
    // second layer
    //

    // Check volume
    ++vol_itr;
    bounds = {64., 80., -500., 500, -M_PI, M_PI};
    range = {602, 1054};
    sf_finder_link = {sf_finder_ids::e_z_phi_grid, 1};

    // Test the links in the volumes
    test_volume_links(vol_itr, 9, bounds, range, sf_finder_link);

    // Check links of modules
    range = {602, 1050};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_rectangle2, 224},
                      {material_ids::e_slab, 602}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1050, 1052};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 30},
                      {material_ids::e_slab, 1050}, portal_mat, {8, 10});

    // disc portals
    range = {1052, 1054};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 24},
                      {material_ids::e_slab, 1052}, portal_mat, {6, 14});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {80., 108., -500., 500, -M_PI, M_PI};
    range = {1054, 1058};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 10, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {1054, 1056};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 32},
                      {material_ids::e_slab, 1054}, portal_mat, {9, 11});
    // disc portals
    range = {1056, 1058};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 26},
                      {material_ids::e_slab, 1056}, portal_mat, {6, 14});

    //
    // third layer
    //

    // Check volume
    ++vol_itr;
    bounds = {108., 124., -500., 500, -M_PI, M_PI};
    range = {1058, 1790};
    sf_finder_link = {sf_finder_ids::e_z_phi_grid, 2};

    // Test the links in the volumes
    test_volume_links(vol_itr, 11, bounds, range, sf_finder_link);

    // Check links of modules
    range = {1058, 1786};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_rectangle2, 672},
                      {material_ids::e_slab, 1058}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1786, 1788};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 34},
                      {material_ids::e_slab, 1786}, portal_mat, {10, 12});

    // disc portals
    range = {1788, 1790};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 28},
                      {material_ids::e_slab, 1788}, portal_mat, {6, 14});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {124., 164., -500., 500, -M_PI, M_PI};
    range = {1790, 1794};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 12, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {1790, 1792};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 36},
                      {material_ids::e_slab, 1790}, portal_mat, {11, 13});
    // disc portals
    range = {1792, 1794};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 30},
                      {material_ids::e_slab, 1792}, portal_mat, {6, 14});

    //
    // fourth layer
    //

    // Check volume
    ++vol_itr;
    bounds = {164., 180., -500., 500, -M_PI, M_PI};
    range = {1794, 2890};
    sf_finder_link = {sf_finder_ids::e_z_phi_grid, 3};

    // Test the links in the volumes
    test_volume_links(vol_itr, 13, bounds, range, sf_finder_link);

    // Check links of modules
    range = {1794, 2886};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_rectangle2, 1400},
                      {material_ids::e_slab, 1794}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {2886, 2888};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 38},
                      {material_ids::e_slab, 2886}, portal_mat,
                      {12, leaving_world});

    // disc portals
    range = {2888, 2890};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 32},
                      {material_ids::e_slab, 2888}, portal_mat, {6, 14});

    //
    // positive endcap
    //

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 500., 595., -M_PI, M_PI};
    range = {2890, 2900};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 14, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {2890, 2892};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 40},
                      {material_ids::e_slab, 2890}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {2892, 2900};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 34},
                      {material_ids::e_slab, 2892}, portal_mat,
                      {15, 7, 8, 9, 10, 11, 12, 13});

    //
    // pos endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 595., 605., -M_PI, M_PI};
    range = {2900, 3012};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 3};

    // Test the links in the volumes
    test_volume_links(vol_itr, 15, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {2900, 3008};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 324},
                      {material_ids::e_slab, 2900}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3008, 3010};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 42},
                      {material_ids::e_slab, 3008}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {3010, 3012};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 42},
                      {material_ids::e_slab, 3010}, portal_mat, {14, 16});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 605., 695., -M_PI, M_PI};
    range = {3012, 3016};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 16, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3012, 3014};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 44},
                      {material_ids::e_slab, 3012}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {3014, 3016};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 44},
                      {material_ids::e_slab, 3014}, portal_mat, {15, 17});

    //
    // pos endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 695., 705., -M_PI, M_PI};
    range = {3016, 3128};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 4};

    // Test the links in the volumes
    test_volume_links(vol_itr, 17, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {3016, 3124};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 432},
                      {material_ids::e_slab, 3016}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3124, 3126};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 46},
                      {material_ids::e_slab, 3124}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {3126, 3128};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 46},
                      {material_ids::e_slab, 3126}, portal_mat, {16, 18});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 705., 815., -M_PI, M_PI};
    range = {3128, 3132};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0};

    // Test the links in the volumes
    test_volume_links(vol_itr, 18, bounds, range, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3128, 3130};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 48},
                      {material_ids::e_slab, 3128}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {3130, 3132};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 48},
                      {material_ids::e_slab, 3130}, portal_mat, {17, 19});

    //
    // pos endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    bounds = {27., 180., 815., 825., -M_PI, M_PI};
    range = {3132, 3244};
    sf_finder_link = {sf_finder_ids::e_r_phi_grid, 5};

    // Test the links in the volumes
    test_volume_links(vol_itr, 19, bounds, range, sf_finder_link);

    // Check the trapezoid modules
    range = {3132, 3240};
    test_module_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_trapezoid2, 540},
                      {material_ids::e_slab, 3132}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3240, 3242};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_cylinder3, 50},
                      {material_ids::e_slab, 3240}, portal_mat,
                      {0, leaving_world});
    // disc portals
    range = {3242, 3244};
    test_portal_links(vol_itr->index(), surfaces.begin() + range[0], range,
                      range[0], {mask_ids::e_portal_ring2, 50},
                      {material_ids::e_slab, 3242}, portal_mat,
                      {18, leaving_world});
}
