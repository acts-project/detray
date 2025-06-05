/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <type_traits>

namespace detray {
/// Call for grid types
template <typename acc_t>
inline void test_finder(const acc_t& finder, const dindex volume_index,
                        const darray<dindex, 2>& range) {

    std::vector<dindex> indices;

    // Check if the correct volume is linked to surface in grid
    // and record the surface indices
    for (const auto& sf : finder.all()) {
        EXPECT_EQ(sf.volume(), volume_index);
        indices.push_back(sf.index());
    }

    // Check if surface indices are consistent with given range of indices
    EXPECT_EQ(finder.size(), range[1] - range[0]);
    for (dindex idx : detray::views::iota(range)) {
        EXPECT_TRUE(std::ranges::find(indices, idx) != indices.end());
    }
}

/// Call for material grid types
template <concepts::scalar scalar_t, typename mat_map_t>
requires std::is_same_v<typename mat_map_t::value_type,
                        material_slab<scalar_t>> inline void
test_mat_map(const mat_map_t& mat_map, const bool is_cyl) {

    // Check if the number of bins is correct
    if (is_cyl) {
        EXPECT_EQ(mat_map.nbins(), 400u);
        EXPECT_EQ(mat_map.size(), 400u);

        auto n_bins = mat_map.axes().nbins_per_axis();
        EXPECT_EQ(n_bins[0], 20u);
        EXPECT_EQ(n_bins[1], 20u);

        auto r_axis = mat_map.template get_axis<axis::label::e_rphi>();
        EXPECT_EQ(r_axis.nbins(), 20u);
        auto z_axis = mat_map.template get_axis<axis::label::e_cyl_z>();
        EXPECT_EQ(z_axis.nbins(), 20u);

        for (const auto& mat_slab : mat_map.all()) {
            EXPECT_TRUE(mat_slab.get_material() ==
                            toy_det_config<scalar_t>{}.mapped_material() ||
                        mat_slab.get_material() == beryllium_tml<scalar_t>{});
        }
    } else {
        EXPECT_EQ(mat_map.nbins(), 60u);
        EXPECT_EQ(mat_map.size(), 60u);

        auto n_bins = mat_map.axes().nbins_per_axis();
        EXPECT_EQ(n_bins[0], 3u);
        EXPECT_EQ(n_bins[1], 20u);

        auto r_axis = mat_map.template get_axis<axis::label::e_r>();
        EXPECT_EQ(r_axis.nbins(), 3u);
        auto z_axis = mat_map.template get_axis<axis::label::e_phi>();
        EXPECT_EQ(z_axis.nbins(), 20u);

        for (const auto& mat_slab : mat_map.all()) {
            EXPECT_TRUE(mat_slab.get_material() ==
                        toy_det_config<scalar_t>{}.mapped_material());
        }
    }
}

template <concepts::algebra algebra_t, typename bfield_t>
inline bool toy_detector_test(
    const detector<toy_metadata<algebra_t>, bfield_t>& toy_det,
    const typename detector<toy_metadata<algebra_t>, bfield_t>::name_map&
        names) {

    using detector_t = detector<toy_metadata<algebra_t>, bfield_t>;
    using scalar_t = dscalar<typename detector_t::algebra_type>;
    using geo_obj_ids = typename detector_t::geo_obj_ids;
    using volume_t = typename detector_t::volume_type;
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using geo_context_t = typename detector_t::geometry_context;
    using mask_ids = typename detector_t::masks::id;
    using mask_link_t = typename detector_t::surface_type::mask_link;
    using material_ids = typename detector_t::materials::id;
    using material_link_t = typename detector_t::surface_type::material_link;
    using accel_ids = typename detector_t::accel::id;
    using accel_link_t = typename volume_t::accel_link_type::index_type;

    EXPECT_EQ(names.at(0u), "toy_detector");

    // Test general consistency
    detail::check_consistency(toy_det, true, names);

    geo_context_t ctx{};
    auto& volumes = toy_det.volumes();
    auto& surfaces = toy_det.surfaces();
    auto& accel = toy_det.accelerator_store();
    auto& transforms = toy_det.transform_store();
    auto& masks = toy_det.mask_store();
    auto& materials = toy_det.material_store();

    // Materials
    auto portal_mat =
        material_slab<scalar_t>(toy_det_config<scalar_t>{}.mapped_material(),
                                1.5f * unit<scalar_t>::mm);
    auto beampipe_mat = material_slab<scalar_t>(beryllium_tml<scalar_t>(),
                                                0.8f * unit<scalar_t>::mm);
    auto pixel_mat = material_slab<scalar_t>(silicon_tml<scalar_t>(),
                                             1.5f * unit<scalar_t>::mm);

    // Link to outer world (leaving detector)
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
    constexpr auto inv_link{detail::invalid_value<dindex>()};
    const bool has_grids =
        (accel.template size<accel_ids::e_cylinder2_grid>() != 0u) ||
        (accel.template size<accel_ids::e_disc_grid>() != 0u);
    const bool has_hom_material =
        (materials.template size<material_ids::e_slab>() != 0);
    const bool has_material_maps =
        (materials.template size<material_ids::e_disc2_map>() != 0);

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 22u);
    EXPECT_EQ(toy_det.surfaces().size(), 3230);
    EXPECT_EQ(transforms.size(ctx), 3252u);
    EXPECT_EQ(masks.template size<mask_ids::e_rectangle2>(), 4u);
    EXPECT_EQ(masks.template size<mask_ids::e_trapezoid2>(), 12u);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_cylinder2>(), 56u);
    EXPECT_EQ(masks.template size<mask_ids::e_portal_ring2>(), 60u);
    EXPECT_EQ(accel.template size<accel_ids::e_brute_force>(), 22u);
    if (has_grids) {
        EXPECT_EQ(accel.template size<accel_ids::e_cylinder2_grid>(), 4);
        EXPECT_EQ(accel.template size<accel_ids::e_disc_grid>(), 6);
    }
    if (has_hom_material) {
        EXPECT_FALSE(has_material_maps);
        EXPECT_EQ(materials.template size<material_ids::e_slab>(), 3141u);
    } else if (has_material_maps) {
        EXPECT_FALSE(has_hom_material);
        EXPECT_EQ(
            materials.template size<material_ids::e_concentric_cylinder2_map>(),
            44u);
        EXPECT_EQ(materials.template size<material_ids::e_disc2_map>(), 42u);
    }

    // Check the surface source links
    /*for (dindex i = 0u; i < toy_det.surfaces().size(); ++i) {

        // The source link was prepared during building as: sf index + 42
        std::uint64_t source{i + 42u};
        const auto& ref_sf = toy_det.surface(i);

        // Only sensitive surfaces were given a source link in the toy detector
        if (ref_sf.is_sensitive()) {
            const auto& result_sf = toy_det.surface(default_searcher{source});

            EXPECT_EQ(result_sf.index(), i) << result_sf;
            EXPECT_EQ(ref_sf, result_sf)
                << "expected: " << ref_sf << ", result: " << result_sf;
        }
    }*/

    /// Test the surface ranges in the volume
    auto check_sf_ranges = [](const typename detector_t::volume_type& vol,
                              dindex_range pt_range, dindex_range sf_range,
                              dindex_range psv_range) {
        EXPECT_EQ(vol.template sf_link<surface_id::e_portal>(), pt_range);
        EXPECT_EQ(vol.template sf_link<surface_id::e_sensitive>(), sf_range);
        EXPECT_EQ(vol.template sf_link<surface_id::e_passive>(), psv_range);
    };

    /// Test the links of a volume.
    ///
    /// @param vol_index volume the modules belong to
    /// @param sf_itr iterator into the surface container, start of the modules
    /// @param range index range of the modules in the surface container
    /// @param trf_index index of the transform (trf container) for the module
    /// @param mask_index type and index of module mask in respective mask cont
    /// @param volume_links links to next volume and next surfaces finder
    auto test_volume_links = [&](decltype(volumes.begin())& vol_itr,
                                 const dindex vol_index,
                                 const darray<dindex, 1>& range,
                                 const accel_link_t& /*accel_link*/) {
        EXPECT_EQ(vol_itr->index(), vol_index);
        EXPECT_EQ(vol_itr->template accel_link<geo_obj_ids::e_portal>().id(),
                  accel_ids::e_brute_force);
        EXPECT_EQ(vol_itr->template accel_link<geo_obj_ids::e_portal>().index(),
                  range[0]);
    };

    /// Test the links of portals (into the next volume or invalid if we leave
    /// the detector).
    ///
    /// @param vol_index volume the portals belong to
    /// @param sf_itr iterator into the surface container, start of the portals
    /// @param range index range of the portals in the surface container
    /// @param trf_index index of the transform (trf container) for the portal
    /// @param mask_index type and index of portal mask in respective mask cont
    /// @param volume_links links to next volume contained in the masks
    auto test_portal_links = [&](const dindex vol_index,
                                 decltype(surfaces.begin())&& sf_itr,
                                 const darray<dindex, 2>& range,
                                 dindex trf_index, mask_link_t&& mask_link,
                                 material_link_t&& material_index,
                                 const material_slab<scalar_t>& mat,
                                 const dvector<dvector<dindex>>&&
                                     volume_links) {
        for (dindex pti = range[0]; pti < range[1]; ++pti) {
            EXPECT_EQ(sf_itr->volume(), vol_index);
            EXPECT_EQ(sf_itr->id(), surface_id::e_portal);
            EXPECT_EQ(sf_itr->index(), pti);
            // The volume index compensates for the number of volume
            // transforms in the transform store
            EXPECT_EQ(sf_itr->transform(), trf_index + vol_index + 1);
            EXPECT_EQ(sf_itr->mask(), mask_link);
            const geometry::surface sf{toy_det, *sf_itr};
            const auto m_volume_links = sf.volume_links();
            EXPECT_EQ(m_volume_links.size(), mask_link.index().size());
            for (std::size_t i = 0u; i < mask_link.index().size(); ++i) {
                EXPECT_EQ(m_volume_links[i], volume_links[pti - range[0]][i])
                    << i;
            }
            if (has_hom_material) {
                EXPECT_EQ(sf_itr->material(), material_index);
                if (sf_itr->material().id() != material_ids::e_none) {
                    EXPECT_EQ(
                        materials.template get<
                            material_ids::e_slab>()[sf_itr->material().index()],
                        mat);
                }
            } else if (has_material_maps) {
                auto mat_link = sf_itr->material();
                if (mat_link.id() == material_ids::e_concentric_cylinder2_map) {
                    test_mat_map<scalar_t>(
                        materials.template get<
                            material_ids::e_concentric_cylinder2_map>()
                            [mat_link.index()],
                        true);
                } else if (mat_link.id() == material_ids::e_disc2_map) {
                    test_mat_map<scalar_t>(
                        materials.template get<
                            material_ids::e_disc2_map>()[mat_link.index()],
                        false);
                }
            } else {
                EXPECT_EQ(sf_itr->material().id(), material_ids::e_none);
            }

            ++sf_itr;
            ++trf_index;
            mask_link.shift(1u);
            if (sf_itr->material().id() != material_ids::e_none) {
                ++material_index;
            }
        }
    };

    /// Test the links of module surface (alway stay in their volume).
    ///
    /// @param vol_index volume the modules belong to
    /// @param sf_itr iterator into the surface container, start of the modules
    /// @param range index range of the modules in the surface container
    /// @param trf_index index of the transform (trf container) for the module
    /// @param mask_index type and index of module mask in respective mask cont
    /// @param volume_links links to next volume contained in the masks
    auto test_module_links =
        [&](const dindex vol_index, decltype(surfaces.begin())&& sf_itr,
            const darray<dindex, 2>& range, dindex trf_index,
            mask_link_t&& mask_index, material_link_t&& material_index,
            const material_slab<scalar_t>& mat,
            const dvector<dindex>&& volume_links,
            bool is_deduplicated = false) {
            for (dindex pti = range[0]; pti < range[1]; ++pti) {
                EXPECT_EQ(sf_itr->volume(), vol_index);
                EXPECT_FALSE(sf_itr->id() == surface_id::e_portal)
                    << sf_itr->barcode();
                EXPECT_EQ(sf_itr->index(), pti);
                // The volume index compensates for the number of volume
                // transforms in the transform store
                EXPECT_EQ(sf_itr->transform(), trf_index + vol_index + 1);
                EXPECT_EQ(sf_itr->mask(), mask_index);
                const geometry::surface sf{toy_det, *sf_itr};
                const auto m_volume_links = sf.volume_links();
                EXPECT_EQ(m_volume_links[0], volume_links[0]);
                if (has_hom_material) {
                    EXPECT_EQ(sf_itr->material(), material_index);
                    EXPECT_EQ(
                        materials.template get<
                            material_ids::e_slab>()[sf_itr->material().index()],
                        mat)
                        << sf_itr->material();
                } else if (has_material_maps &&
                           (sf_itr->id() == surface_id::e_passive)) {
                    // beampipe
                    auto mat_link = sf_itr->material();
                    EXPECT_EQ(mat_link.id(),
                              material_ids::e_concentric_cylinder2_map);

                    test_mat_map<scalar_t>(
                        materials.template get<
                            material_ids::e_concentric_cylinder2_map>()
                            [mat_link.index()],
                        true);
                } else {
                    EXPECT_EQ(sf_itr->material().id(), material_ids::e_none);
                }

                ++sf_itr;
                ++trf_index;
                if (!is_deduplicated) {
                    ++mask_index;
                }
                ++material_index;
            }
        };

    /// Test the detectors acceleration data structures.
    ///
    /// @param vol_itr iterator over the volume descriptors
    /// @param accel_store the detectors acceleration data structure store
    /// @param pt_range index range of the portals in the surface lookup
    /// @param sf_range index range of the surfaces in the surface lookup
    auto test_accel =
        [has_grids](
            decltype(volumes.begin())& vol_itr,
            const typename detector_t::accelerator_container& accel_store,
            const darray<dindex, 2>& pt_range,
            const darray<dindex, 2>& sf_range = {0u, 0u}) {
            // Link to the acceleration data structures the volume holds
            const auto& link = vol_itr->accel_link();

            // Test the portal search
            const auto& bf_finder =
                accel_store
                    .template get<accel_ids::e_brute_force>()[link[0].index()];

            // This means no grids, all surfaces are in the brute force method
            if (!has_grids) {
                const auto full_range = darray<dindex, 2>{
                    sf_range[0], math::max(pt_range[1], sf_range[1])};
                test_finder(bf_finder, vol_itr->index(), full_range);
            } else {
                test_finder(bf_finder, vol_itr->index(), pt_range);

                // Test the module search if grids were filled
                if (!link[1].is_invalid()) {
                    if (link[1].id() == accel_ids::e_cylinder2_grid) {
                        const auto& cyl_grid = accel_store.template get<
                            accel_ids::e_cylinder2_grid>()[link[1].index()];
                        test_finder(cyl_grid, vol_itr->index(), sf_range);
                    } else {
                        const auto& disc_grid = accel_store.template get<
                            accel_ids::e_disc_grid>()[link[1].index()];
                        test_finder(disc_grid, vol_itr->index(), sf_range);
                    }
                }
            }
        };

    //
    // beampipe
    //

    // Check volume
    auto vol_itr = volumes.begin();
    EXPECT_EQ(names.at(vol_itr->index() + 1), "beampipe_0");
    darray<dindex, 1> index = {0u};
    accel_link_t accel_link{accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {1u, 6u}, {}, {0u, 1u});

    // Test the links in the volumes
    test_volume_links(vol_itr, 0u, index, accel_link);

    // Check links of beampipe itself
    darray<dindex, 2> range = {0u, 1u};
    test_module_links(vol_itr->index(), surfaces.begin(), range, range[0],
                      {mask_ids::e_cylinder2, {0u, 1u}},
                      {material_ids::e_slab, 0u}, beampipe_mat,
                      {vol_itr->index()});

    // Check links of portals
    // left disc portal
    range = {1u, 2u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {0u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{leaving_world}});
    // cylinder portals (neg. endcap)
    range = {2u, 3u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, {1u, 6u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{1u, 2u, 3u, 4u, 5u, 6u}});
    // central cylinder portal
    range = {3u, 4u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, {7u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat, {{8u}});
    // right disc portal
    range = {4u, 5u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {1u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{leaving_world}});
    // cylinder portals (pos. endcap)
    range = {5u, 6u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_cylinder2, {8u, 6u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{16u, 17u, 18u, 19u, 20u, 21u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {0u, 6u});

    //
    // neg endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_1");
    range = {6u, 118u};
    index = {1u};
    accel_link = {accel_ids::e_disc_grid, 0u};
    check_sf_ranges(*vol_itr, {114u, 118u}, {6u, 114u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 1u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {6u, 46u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {0u, 1u}},
                      {material_ids::e_slab, 1u}, pixel_mat, {vol_itr->index()},
                      true);
    // One mask for the outer ring
    range = {46u, 114u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {1u, 1u}},
                      {material_ids::e_slab, 41u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {114u, 116u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {14u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {116u, 118u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {2u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{4u}, {2u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {114u, 118u}, {6u, 114u});

    //
    // connector gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "connector_gap_2");
    range = {118u, 122u};
    index = {2u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {118u, 122u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 2u, index, accel_link);

    // Check links of portals
    // disk portals
    range = {118u, 119u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {4u, 9u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u}});
    range = {119u, 120u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {13u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat, {{1u}});
    // cylinder portals
    range = {120u, 122u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {16u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {118u, 122u});

    //
    // neg endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_3");
    range = {122u, 234u};
    index = {3u};
    accel_link = {accel_ids::e_disc_grid, 1u};
    check_sf_ranges(*vol_itr, {230u, 234u}, {122u, 230u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 3u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {122u, 162u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {2u, 1u}},
                      {material_ids::e_slab, 109u}, pixel_mat,
                      {vol_itr->index()}, true);
    // One mask for the outer ring
    range = {162u, 230u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {3u, 1u}},
                      {material_ids::e_slab, 149u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {230u, 232u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {18u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {232u, 234u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {14u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{6u}, {4u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {230u, 234u}, {122u, 230u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_4");
    range = {234u, 238u};
    index = {4u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {234u, 238u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 4u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {234u, 236u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {20u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {236u, 238u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {16u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{3u}, {1u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {234u, 238u});

    //
    // neg endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_5");
    range = {238u, 350u};
    index = {5u};
    accel_link = {accel_ids::e_disc_grid, 2u};
    check_sf_ranges(*vol_itr, {346u, 350u}, {238u, 346u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 5u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {238u, 278u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {4u, 1u}},
                      {material_ids::e_slab, 217u}, pixel_mat,
                      {vol_itr->index()}, true);
    // One mask for the outer ring
    range = {278u, 346u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {5u, 1u}},
                      {material_ids::e_slab, 257u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {346u, 348u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {22u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {348u, 350u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {18u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{leaving_world}, {6u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {346u, 350u}, {238u, 346u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_6");
    range = {350u, 354u};
    index = {6u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {350u, 354u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 6u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {350u, 352u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {24u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {352u, 354u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {20u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{5u}, {3u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {350u, 354u});

    //
    // barrel
    //

    //
    // first layer
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "barrel_7");
    range = {354u, 582u};
    index = {7u};
    accel_link = {accel_ids::e_cylinder2_grid, 0u};
    check_sf_ranges(*vol_itr, {578u, 582u}, {354u, 578u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 7u, index, accel_link);

    // Check links of modules
    range = {354u, 578u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, {0u, 1u}},
                      {material_ids::e_slab, 325u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {578u, 580u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {26u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{8u}, {10u}});

    // disc portals
    range = {580u, 582u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {22u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {578u, 582u}, {354u, 578u});

    //
    // gap (between first barrel layer and beampipe volume)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_8");
    range = {582u, 584u};
    index = {8u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {582u, 586u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 8u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {582u, 584u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {28u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {7u}});
    // disc portals
    range = {584u, 586u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {24u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {582u, 586u});

    //
    // second layer
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "barrel_9");
    range = {586u, 1038u};
    index = {9u};
    accel_link = {accel_ids::e_cylinder2_grid, 1u};
    check_sf_ranges(*vol_itr, {1034u, 1038u}, {586u, 1034u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 9u, index, accel_link);

    // Check links of modules
    range = {586u, 1034u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, {1u, 1u}},
                      {material_ids::e_slab, 549u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {1034u, 1036u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {30u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{10u}, {12u}});

    // disc portals
    range = {1036u, 1038u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {26u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {1034u, 1038u}, {586u, 1034u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_10");
    range = {1038u, 1042u};
    index = {10u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {1038u, 1042u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 10u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {1038u, 1040u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {32u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{7u}, {9u}});
    // disc portals
    range = {1040u, 1042u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {28u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {1038u, 1042u});

    //
    // third layer
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "barrel_11");
    range = {1042u, 1774u};
    index = {11u};
    accel_link = {accel_ids::e_cylinder2_grid, 2u};
    check_sf_ranges(*vol_itr, {1770u, 1774u}, {1042u, 1770u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 11u, index, accel_link);

    // Check links of modules
    range = {1042u, 1770u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, {2u, 1u}},
                      {material_ids::e_slab, 997u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {1770u, 1772u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {34u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{12u}, {14u}});

    // disc portals
    range = {1772u, 1774u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {30u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {1770u, 1774u}, {1042u, 1770u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_12");
    range = {1774u, 1778u};
    index = {12u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {1774u, 1778u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 12u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {1774u, 1776u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {36u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{9u}, {11u}});
    // disc portals
    range = {1776u, 1778u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {32u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {1774u, 1778u});

    //
    // fourth layer
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "barrel_13");
    range = {1778u, 2874u};
    index = {13u};
    accel_link = {accel_ids::e_cylinder2_grid, 3u};
    check_sf_ranges(*vol_itr, {2870u, 2874u}, {1778u, 2870u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 13u, index, accel_link);

    // Check links of modules
    range = {1778u, 2870u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_rectangle2, {3u, 1u}},
                      {material_ids::e_slab, 1725u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {2870u, 2872u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {38u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{14u}, {15u}});

    // disc portals
    range = {2872u, 2874u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {34u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {2870u, 2874u}, {1778u, 2870u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_14");
    range = {2874u, 2878u};
    index = {14u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {2874u, 2878u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 14u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {2874u, 2876u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {40u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{11u}, {13u}});
    // disc portals
    range = {2876u, 2878u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {36u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {2874u, 2878u});

    //
    // gap (between last barrel layer and full detector outer radius)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_15");
    range = {2878u, 2882u};
    index = {15u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {2878u, 2882u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 15u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {2878u, 2880u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {42u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{13u}, {leaving_world}});
    // disc portals
    range = {2880u, 2882u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {38u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{2u}, {17u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {2878u, 2882u});

    //
    // positive endcap
    //

    //
    // pos endcap (layer 1)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_16");
    range = {2882u, 2994u};
    index = {16u};
    accel_link = {accel_ids::e_disc_grid, 0u};
    check_sf_ranges(*vol_itr, {2990u, 2994u}, {2882u, 2990u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 16u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {2882u, 2922u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {6u, 1u}},
                      {material_ids::e_slab, 2817u}, pixel_mat,
                      {vol_itr->index()}, true);
    // One mask for the outer ring
    range = {2922u, 2990u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {7u, 1u}},
                      {material_ids::e_slab, 2857u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {2990u, 2992u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {44u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {2992u, 2994u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {40u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{17u}, {19u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {2990u, 2994u}, {2882u, 2990u});

    //
    // connector gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "connector_gap_17");
    range = {2994u, 2998u};
    index = {17u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {2994u, 2998u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 17u, index, accel_link);

    // Check links of portals
    // disk portals
    range = {2994u, 2995u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {42u, 9u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{7u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u}});
    range = {2995u, 2996u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {51u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat, {{16u}});
    // cylinder portals
    range = {2996u, 2998u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {46u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {2994u, 2998u});

    //
    // pos endcap (layer 2)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_18");
    range = {2998u, 3110u};
    index = {18u};
    accel_link = {accel_ids::e_disc_grid, 3u};
    check_sf_ranges(*vol_itr, {3106u, 3110u}, {2998u, 3106u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 18u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {2998u, 3038u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {8u, 1u}},
                      {material_ids::e_slab, 2925u}, pixel_mat,
                      {vol_itr->index()}, true);
    // One mask for the outer ring
    range = {3038u, 3106u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {9u, 1u}},
                      {material_ids::e_slab, 2965u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {3106u, 3108u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {48u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {3108u, 3110u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {52u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{19u}, {21u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {3106u, 3110u}, {2998u, 3106u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_19");
    range = {3110u, 3114u};
    index = {19u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {3110u, 3114u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 19u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {3110u, 3112u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {50u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {3112u, 3114u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {54u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{16u}, {18u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {3110u, 3114u});

    //
    // pos endcap (layer 3)
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "endcap_20");
    range = {3114u, 3226u};
    index = {20u};
    accel_link = {accel_ids::e_disc_grid, 4u};
    check_sf_ranges(*vol_itr, {3222u, 3226u}, {3114u, 3222u}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 20u, index, accel_link);

    // Check the trapezoid modules
    // One mask for the inner ring
    range = {3114u, 3154u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {10u, 1u}},
                      {material_ids::e_slab, 3033u}, pixel_mat,
                      {vol_itr->index()}, true);
    // One mask for the outer ring
    range = {3154u, 3222u};
    test_module_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_trapezoid2, {11u, 1u}},
                      {material_ids::e_slab, 3073u}, pixel_mat,
                      {vol_itr->index()}, true);

    // Check links of portals
    // cylinder portals
    range = {3222u, 3224u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {52u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {3224u, 3226u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {56u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{21u}, {leaving_world}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {3222u, 3226u}, {3114u, 3222u});

    //
    // gap
    //

    // Check volume
    ++vol_itr;
    EXPECT_EQ(names.at(vol_itr->index() + 1), "gap_21");
    range = {3226u, 3230u};
    index = {21u};
    accel_link = {accel_ids::e_brute_force, 0u};
    check_sf_ranges(*vol_itr, {3226u, 3230u}, {}, {});

    // Test the links in the volumes
    test_volume_links(vol_itr, 21u, index, accel_link);

    // Check links of portals
    // cylinder portals
    range = {3226u, 3228u};
    test_portal_links(
        vol_itr->index(),
        surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]), range,
        range[0], {mask_ids::e_portal_cylinder2, {54u, 1u}},
        {material_ids::e_none, inv_link}, portal_mat, {{0u}, {leaving_world}});
    // disc portals
    range = {3228u, 3230u};
    test_portal_links(vol_itr->index(),
                      surfaces.begin() + static_cast<std::ptrdiff_t>(range[0]),
                      range, range[0], {mask_ids::e_portal_ring2, {58u, 1u}},
                      {material_ids::e_none, inv_link}, portal_mat,
                      {{18u}, {20u}});

    // Check link of surfaces in surface finder
    test_accel(vol_itr, accel, {3226u, 3230u});

    return true;
}

}  // namespace detray
