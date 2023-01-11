/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#include <gtest/gtest.h>

#include <type_traits>

#include "detray/detectors/create_toy_geometry.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace detray;

namespace {

/// Functor that returns the volume/sf_finder links of a particluar mask
/// instance in the mask container of a detector
struct volume_link_getter {

    template <typename mask_group_t, typename mask_range_t>
    inline auto operator()(mask_group_t& mask_group, mask_range_t& mask_range) {
        return mask_group[mask_range].volume_link();
    }
};

/// Functor that tests the volume indices of surfaces in the grid
struct surface_grid_tester {

    /// Call for grid types
    template <typename grid_group_t, typename surface_container_t,
              std::enable_if_t<
                  not std::is_same_v<typename grid_group_t::grid_type, void>,
                  bool> = true>
    inline void operator()(const grid_group_t& grid_group,
                           const dindex grid_index, const dindex volume_index,
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
    }

    /// Call for the brute force type (do nothing)
    template <
        typename grid_group_t, typename surface_container_t,
        std::enable_if_t<std::is_same_v<typename grid_group_t::grid_type, void>,
                         bool> = true>
    inline bool operator()(const grid_group_t&, const dindex, const dindex,
                           const surface_container_t&,
                           const darray<dindex, 2>&) {
        return true;
    }
};

}  // anonymous namespace

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    vecmem::host_memory_resource host_mr;
    constexpr std::size_t n_brl_layers{4u};
    constexpr std::size_t n_edc_layers{3u};

    const auto toy_det =
        create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    using detector_t = decltype(toy_det);
    using geo_obj_ids = typename detector_t::geo_obj_ids;
    using volume_t = typename detector_t::volume_type;
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using geo_context_t = typename detector_t::geometry_context;
    using mask_ids = typename detector_t::masks::id;
    using mask_link_t = typename detector_t::surface_type::mask_link;
    using material_ids = typename detector_t::materials::id;
    using material_link_t = typename detector_t::surface_type::material_link;
    using sf_finder_ids = typename detector_t::sf_finders::id;
    using sf_finder_link_t = typename volume_t::link_type::index_type;

    geo_context_t ctx{};
    auto& volumes = toy_det.volumes();
    auto& surfaces = toy_det.surfaces();
    // auto& sf_finders = toy_det.sf_finder_store();
    auto& transforms = toy_det.transform_store();
    auto& masks = toy_det.mask_store();
    auto& materials = toy_det.material_store();

    // Materials
    auto portal_mat =
        material_slab<scalar>(vacuum<scalar>(), 0.f * unit<scalar>::mm);
    auto beampipe_mat =
        material_slab<scalar>(beryllium_tml<scalar>(), 0.8f * unit<scalar>::mm);
    auto pixel_mat =
        material_slab<scalar>(silicon_tml<scalar>(), 0.15f * unit<scalar>::mm);

    // Link to outer world (leaving detector)
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 20u);
    EXPECT_EQ(surfaces.size(), 3244u);
    /*EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_brute_force>(), 1);
    EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_cylinder_grid>(),
              n_brl_layers);
    EXPECT_EQ(sf_finders.template size<sf_finder_ids::e_disc_grid>(), 14);*/
    EXPECT_EQ(transforms.size(ctx), 3244u);
    EXPECT_EQ(masks.template size<mask_ids::e_rectangle2>(), 2492u);
    EXPECT_EQ(masks.template size<mask_ids::e_trapezoid2>(), 648u);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_cylinder2>(), 52u);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_ring2>(), 52u);
    EXPECT_EQ(materials.template size<material_ids::e_slab>(), 3244u);

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
            const darray<scalar, 6>& bounds, const darray<dindex, 1>& range,
            const sf_finder_link_t& /*sf_finder_link*/) {
            EXPECT_EQ(vol_itr->index(), vol_index);
            EXPECT_EQ(vol_itr->bounds(), bounds);
            EXPECT_EQ(vol_itr->template link<geo_obj_ids::e_portal>().id(),
                      sf_finder_ids::e_brute_force);
            EXPECT_EQ(vol_itr->template link<geo_obj_ids::e_portal>().index(),
                      range[0]);
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
            EXPECT_EQ(sf_itr->id(), surface_id::e_portal);
            EXPECT_EQ(sf_itr->index(), pti);
            EXPECT_EQ(sf_itr->transform(), trf_index);
            EXPECT_EQ(sf_itr->mask(), mask_link);
            const auto volume_link =
                masks.template visit<volume_link_getter>(sf_itr->mask());
            EXPECT_EQ(volume_link, volume_links[pti - range[0]]);
            EXPECT_EQ(
                materials
                    .get<material_ids::e_slab>()[sf_itr->material().index()],
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
            EXPECT_FALSE(sf_itr->id() == surface_id::e_portal);
            EXPECT_EQ(sf_itr->index(), pti);
            EXPECT_EQ(sf_itr->transform(), trf_index);
            EXPECT_EQ(sf_itr->mask(), mask_index);
            EXPECT_EQ(sf_itr->material(), material_index);
            const auto volume_link =
                masks.template visit<volume_link_getter>(sf_itr->mask());
            EXPECT_EQ(volume_link, volume_links[0]);
            EXPECT_EQ(
                materials
                    .get<material_ids::e_slab>()[sf_itr->material().index()],
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
            sf_finder_cont.template visit<surface_grid_tester>(
                vol_itr->sf_finder_link(), vol_itr->index(), surface_cont,
                range);
        };*/

    //
    // beampipe
    //

    // Check volume
    auto vol_itr = volumes.begin();
    darray<scalar, 6> bounds = {
        0.f, 27.f, -825.f, 825.f, -constant<scalar>::pi, constant<scalar>::pi};
    darray<dindex, 1> index = {0u};
    sf_finder_link_t sf_finder_link{sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 0u, bounds, index, sf_finder_link);

    // Check links of beampipe itself
    darray<dindex, 2> range = {0u, 1u};
    test_module_links(vol_itr->index(), surfaces.begin(), range, range[0],
                      {mask_ids::e_cylinder2, 0u}, {material_ids::e_slab, 0u},
                      beampipe_mat, {vol_itr->index()});

    // Check links of portals
    // cylinder portals
    range = {1u, 8u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 1u},
                      {material_ids::e_slab, 1u}, portal_mat,
                      {1u, 2u, 3u, 4u, 5u, 6u, 7u});
    range = {8u, 14u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 8u},
                      {material_ids::e_slab, 8u}, portal_mat,
                      {14u, 15u, 16u, 17u, 18u, 19u});

    // disc portals
    range = {14u, 16u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 0u},
                      {material_ids::e_slab, 14u}, portal_mat,
                      {leaving_world, leaving_world});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -825.f,
              -815.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {16u, 128u};
    index = {1u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 1u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {16u, 124u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 0u},
                      {material_ids::e_slab, 16u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {124u, 126u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 14u},
                      {material_ids::e_slab, 124u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {126u, 128u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 2u},
                      {material_ids::e_slab, 126u}, portal_mat,
                      {leaving_world, 2u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -815.f,
              -705.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {128u, 132u};
    index = {2u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 2u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {128u, 130u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 16u},
                      {material_ids::e_slab, 128u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {130u, 132u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 4u},
                      {material_ids::e_slab, 130u}, portal_mat, {1u, 3u});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -705.f,
              -695.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {132u, 244u};
    index = {3u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 1u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 3u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {132u, 240u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 108u},
                      {material_ids::e_slab, 132u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {240u, 242u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 18u},
                      {material_ids::e_slab, 240u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {242u, 244u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 6u},
                      {material_ids::e_slab, 242u}, portal_mat, {2u, 4u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -695.f,
              -605.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {244u, 248u};
    index = {4u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 4u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {244u, 246u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 20u},
                      {material_ids::e_slab, 244u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {246u, 248u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 8u},
                      {material_ids::e_slab, 246u}, portal_mat, {3u, 5u});

    //
    // neg endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -605.f,
              -595.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {248u, 360u};
    index = {5u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 2u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 5u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {248u, 356u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 216u},
                      {material_ids::e_slab, 248u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {356u, 358u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 22u},
                      {material_ids::e_slab, 356u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {358u, 360u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 10u},
                      {material_ids::e_slab, 358u}, portal_mat, {4u, 6u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {27.f,
              180.f,
              -595.f,
              -500.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {360u, 370u};
    index = {6u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 6u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {360u, 362u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 24u},
                      {material_ids::e_slab, 360u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {362u, 370u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 12u},
                      {material_ids::e_slab, 362u}, portal_mat,
                      {5u, 7u, 8u, 9u, 10u, 11u, 12u, 13u});

    //
    // barrel
    //

    //
    // first layer
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 38.f, -500.f, 500.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {370u, 598u};
    index = {7u};
    sf_finder_link = {sf_finder_ids::e_cylinder_grid, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 7u, bounds, index, sf_finder_link);

    // Check links of modules
    range = {370u, 594u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, 0u},
                      {material_ids::e_slab, 370u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {594u, 596u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 26u},
                      {material_ids::e_slab, 594u}, portal_mat, {0u, 8u});

    // disc portals
    range = {596u, 598u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 20u},
                      {material_ids::e_slab, 596u}, portal_mat, {6u, 14u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {
        38.f, 64.f, -500.f, 500.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {598u, 602u};
    index = {8u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 8u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {598u, 600u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 28u},
                      {material_ids::e_slab, 598u}, portal_mat, {7u, 9u});
    // disc portals
    range = {600u, 602u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 22u},
                      {material_ids::e_slab, 600u}, portal_mat, {6u, 14u});

    //
    // second layer
    //

    // Check volume
    ++vol_itr;
    bounds = {
        64.f, 80.f, -500.f, 500.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {602u, 1054u};
    index = {9u};
    sf_finder_link = {sf_finder_ids::e_cylinder_grid, 1u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 9u, bounds, index, sf_finder_link);

    // Check links of modules
    range = {602u, 1050u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, 224u},
                      {material_ids::e_slab, 602u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1050u, 1052u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 30u},
                      {material_ids::e_slab, 1050u}, portal_mat, {8u, 10u});

    // disc portals
    range = {1052u, 1054u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 24u},
                      {material_ids::e_slab, 1052u}, portal_mat, {6u, 14u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {80.f,
              108.f,
              -500.f,
              500.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {1054u, 1058u};
    index = {10u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 10u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {1054u, 1056u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 32u},
                      {material_ids::e_slab, 1054u}, portal_mat, {9u, 11u});
    // disc portals
    range = {1056u, 1058u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 26u},
                      {material_ids::e_slab, 1056u}, portal_mat, {6u, 14u});

    //
    // third layer
    //

    // Check volume
    ++vol_itr;
    bounds = {108.f,
              124.f,
              -500.f,
              500.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {1058u, 1790u};
    index = {11u};
    sf_finder_link = {sf_finder_ids::e_cylinder_grid, 2u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 11u, bounds, index, sf_finder_link);

    // Check links of modules
    range = {1058u, 1786u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, 672u},
                      {material_ids::e_slab, 1058u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {1786u, 1788u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 34u},
                      {material_ids::e_slab, 1786u}, portal_mat, {10u, 12u});

    // disc portals
    range = {1788u, 1790u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 28u},
                      {material_ids::e_slab, 1788u}, portal_mat, {6u, 14u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {124.f,
              164.f,
              -500.f,
              500.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {1790u, 1794u};
    index = {12u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 12u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {1790u, 1792u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 36u},
                      {material_ids::e_slab, 1790u}, portal_mat, {11u, 13u});
    // disc portals
    range = {1792u, 1794u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 30u},
                      {material_ids::e_slab, 1792u}, portal_mat, {6u, 14u});

    //
    // fourth layer
    //

    // Check volume
    ++vol_itr;
    bounds = {164.f,
              180.f,
              -500.f,
              500.f,
              -constant<scalar>::pi,
              constant<scalar>::pi};
    range = {1794u, 2890u};
    index = {13u};
    sf_finder_link = {sf_finder_ids::e_cylinder_grid, 3u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 13u, bounds, index, sf_finder_link);

    // Check links of modules
    range = {1794u, 2886u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, 1400u},
                      {material_ids::e_slab, 1794u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {2886u, 2888u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 38u},
                      {material_ids::e_slab, 2886u}, portal_mat,
                      {12u, leaving_world});

    // disc portals
    range = {2888u, 2890u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 32u},
                      {material_ids::e_slab, 2888u}, portal_mat, {6u, 14u});

    //
    // positive endcap
    //

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 500.f, 595.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {2890u, 2900u};
    index = {14u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 14u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {2890u, 2892u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 40u},
                      {material_ids::e_slab, 2890u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {2892u, 2900u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 34u},
                      {material_ids::e_slab, 2892u}, portal_mat,
                      {15u, 7u, 8u, 9u, 10u, 11u, 12u, 13u});

    //
    // pos endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 595.f, 605.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {2900u, 3012u};
    index = {15u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 3u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 15u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {2900u, 3008u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 324u},
                      {material_ids::e_slab, 2900u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3008u, 3010u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 42u},
                      {material_ids::e_slab, 3008u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {3010u, 3012u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 42u},
                      {material_ids::e_slab, 3010u}, portal_mat, {14u, 16u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 605.f, 695.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {3012u, 3016u};
    index = {16u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 16u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3012u, 3014u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 44u},
                      {material_ids::e_slab, 3012u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {3014u, 3016u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 44u},
                      {material_ids::e_slab, 3014u}, portal_mat, {15u, 17u});

    //
    // pos endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 695.f, 705.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {3016u, 3128u};
    index = {17u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 4u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 17u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {3016u, 3124u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 432u},
                      {material_ids::e_slab, 3016u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3124u, 3126u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 46u},
                      {material_ids::e_slab, 3124u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {3126u, 3128u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 46u},
                      {material_ids::e_slab, 3126u}, portal_mat, {16u, 18u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 705.f, 815.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {3128u, 3132u};
    index = {18u};
    sf_finder_link = {sf_finder_ids::e_brute_force, 0u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 18u, bounds, index, sf_finder_link);

    // Check links of portals
    // cylinder portals
    range = {3128u, 3130u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 48u},
                      {material_ids::e_slab, 3128u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {3130u, 3132u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 48u},
                      {material_ids::e_slab, 3130u}, portal_mat, {17u, 19u});

    //
    // pos endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    bounds = {
        27.f, 180.f, 815.f, 825.f, -constant<scalar>::pi, constant<scalar>::pi};
    range = {3132u, 3244u};
    index = {19u};
    sf_finder_link = {sf_finder_ids::e_disc_grid, 5u};

    // Test the links in the volumes
    test_volume_links(vol_itr, 19u, bounds, index, sf_finder_link);

    // Check the trapezoid modules
    range = {3132u, 3240u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, 540u},
                      {material_ids::e_slab, 3132u}, pixel_mat,
                      {vol_itr->index()});

    // Check link of surfaces in surface finder
    // test_surfaces_grid(vol_itr, sf_finders, surfaces, range);

    // Check links of portals
    // cylinder portals
    range = {3240u, 3242u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, 50u},
                      {material_ids::e_slab, 3240u}, portal_mat,
                      {0u, leaving_world});
    // disc portals
    range = {3242u, 3244u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, 50u},
                      {material_ids::e_slab, 3242u}, portal_mat,
                      {18u, leaving_world});
}
