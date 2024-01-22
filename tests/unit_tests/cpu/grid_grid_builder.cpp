/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/test/types.hpp"
#include "detray/tools/grid_builder.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Gtest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;
using namespace detray::n_axis;

namespace {

using point3 = test::point3;
using vector3 = test::vector3;

test::transform3 Identity{};

using detector_t = detector<toy_metadata>;

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_tools, grid_factory) {

    // Data-owning grid collection
    vecmem::host_memory_resource host_mr;
    auto gr_factory =
        grid_factory<bins::static_array<dindex, 3>, simple_serializer>{host_mr};

    // Build from existing mask of a surface
    const scalar minR{0.f};
    const scalar maxR{10.f};
    const scalar minPhi{0.f};
    const scalar maxPhi{constant<scalar>::pi};
    mask<annulus2D<>> ann2{0u, minR, maxR, minPhi, maxPhi, 0.f, 0.f, 0.f};

    // Grid with correctly initialized axes, but empty bin content
    auto ann_gr = gr_factory.new_grid(ann2, {5, 10});

    // Test axis
    const auto& ann_axis_r = ann_gr.template get_axis<label::e_r>();
    EXPECT_EQ(ann_axis_r.label(), label::e_r);
    EXPECT_EQ(ann_axis_r.bounds(), bounds::e_closed);
    EXPECT_EQ(ann_axis_r.binning(), binning::e_regular);
    EXPECT_EQ(ann_axis_r.nbins(), 5u);
    EXPECT_NEAR(ann_axis_r.span()[0], 0.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(ann_axis_r.span()[1], 10.f,
                std::numeric_limits<scalar>::epsilon());

    // Test fill a bin to see, if bin content was correctly initialized
    point3 p = {0.5f, 2.f, 0.f};
    vector3 d{};
    auto loc_p = ann_gr.project(Identity, p, d);
    ann_gr.template populate<attach<>>(loc_p, 3u);
    ann_gr.template populate<attach<>>(loc_p, 5u);
    auto bin2 = ann_gr.search(loc_p);

    EXPECT_TRUE(bin2.size() == 2u);
    EXPECT_FALSE(bin2.empty());
    EXPECT_EQ(bin2[0], 3u);
    EXPECT_EQ(bin2[1], 5u);

    // gr_builder.to_string(ann_gr);

    // Build from parameters
    const std::vector<scalar> bin_edges_z{-10.f, -8.f, -6.5f, -1.f,
                                          4.f,   5.f,  6.f,   9.f};
    const std::vector<scalar> bin_edges_phi{};

    auto cyl_gr = gr_factory.template new_grid<cylinder2D<>>(
        {bin_edges_z.front(), bin_edges_z.back(), 0.f,
         2.f * constant<scalar>::pi},
        {bin_edges_z.size() - 1, 10u}, {bin_edges_phi, bin_edges_z},
        types::list<circular<label::e_rphi>, closed<label::e_cyl_z>>{},
        types::list<regular<>, irregular<>>{});

    // Test axis
    const auto& cyl_axis_z = cyl_gr.template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
    EXPECT_EQ(cyl_axis_z.binning(), binning::e_irregular);
    EXPECT_EQ(cyl_axis_z.nbins(), 7u);
    EXPECT_NEAR(cyl_axis_z.span()[0], -10.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(cyl_axis_z.span()[1], 9.f,
                std::numeric_limits<scalar>::epsilon());

    // Test fill a bin to see, if bin content was correctly initialized
    loc_p = cyl_gr.project(Identity, p, d);
    cyl_gr.template populate<attach<>>(loc_p, 33u);
    cyl_gr.template populate<attach<>>(loc_p, 55u);

    EXPECT_EQ(cyl_gr.search(loc_p)[0], 33u);
    EXPECT_EQ(cyl_gr.search(loc_p)[1], 55u);

    // Build the same cylinder grid from a mask
    const scalar r{5.f};
    const scalar n_half_z{-10.f};
    const scalar p_half_z{9.f};
    mask<cylinder2D<>> cyl2{0u, r, n_half_z, p_half_z};

    auto cyl_gr2 = gr_factory.template new_grid<circular<label::e_rphi>,
                                                closed<label::e_cyl_z>,
                                                regular<>, irregular<>>(
        cyl2, {bin_edges_z.size() - 1, 10u}, {bin_edges_phi, bin_edges_z});
}

/// Unittest: Test the grid builder
GTEST_TEST(detray_tools, grid_builder) {

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             bins::static_array<detector_t::surface_type, 1>>;

    auto gbuilder =
        grid_builder<detector_t, cyl_grid_t, detray::detail::bin_associator>{
            nullptr};

    // The cylinder portals are at the end of the surface range by construction
    const auto cyl_mask = mask<cylinder2D<>>{0u, 10.f, -500.f, 500.f};
    std::size_t n_phi_bins{5u}, n_z_bins{4u};

    // Build empty grid
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});
    auto cyl_grid = gbuilder.get();

    const auto& cyl_axis_z = cyl_grid.template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
    EXPECT_EQ(cyl_axis_z.binning(), binning::e_regular);
    EXPECT_EQ(cyl_axis_z.nbins(), 4u);
    EXPECT_NEAR(cyl_axis_z.span()[0], -500.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(cyl_axis_z.span()[1], 500.f,
                std::numeric_limits<scalar>::epsilon());
}

/// Integration test: grid builder as volume builder decorator
GTEST_TEST(detray_tools, decorator_grid_builder) {

    using transform3 = typename detector_t::transform3;
    using geo_obj_id = typename detector_t::geo_obj_ids;
    using acc_ids = typename detector_t::accel::id;
    using mask_id = typename detector_t::masks::id;

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             bins::static_array<detector_t::surface_type, 1>>;

    using pt_cylinder_t = cylinder2D<false, cylinder_portal_intersector>;
    using pt_cylinder_factory_t = surface_factory<detector_t, pt_cylinder_t>;
    using rectangle_factory = surface_factory<detector_t, rectangle2D<>>;
    using trapezoid_factory = surface_factory<detector_t, trapezoid2D<>>;
    using cylinder_factory = surface_factory<detector_t, pt_cylinder_t>;

    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};
    const auto vol_idx{
        static_cast<typename detector_t::surface_type::navigation_link>(
            d.volumes().size())};

    auto vbuilder = std::make_unique<volume_builder<detector_t>>(
        volume_id::e_cylinder, vol_idx);
    auto gbuilder = grid_builder<detector_t, cyl_grid_t>{std::move(vbuilder)};
    // passive surfaces are added to the grid
    // gbuilder.set_add_surfaces();

    // The cylinder portals are at the end of the surface range by construction
    const auto cyl_mask = mask<cylinder2D<>>{0u, 10.f, -500.f, 500.f};
    std::size_t n_phi_bins{5u}, n_z_bins{4u};

    // Build empty grid
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});

    EXPECT_TRUE(d.volumes().size() == 0);

    // Add some portals first
    auto pt_cyl_factory = std::make_shared<pt_cylinder_factory_t>();

    typename pt_cylinder_factory_t::sf_data_collection cyl_sf_data;
    cyl_sf_data.emplace_back(surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 0u,
                             std::vector<scalar>{10.f, -1500.f, 1500.f});
    cyl_sf_data.emplace_back(surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1500.f});
    pt_cyl_factory->push_back(std::move(cyl_sf_data));

    // Then some passive and sensitive surfaces
    auto rect_factory = std::make_shared<rectangle_factory>();

    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.emplace_back(surface_id::e_sensitive,
                              transform3(point3{7.07f, 7.07f, -500.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(surface_id::e_sensitive,
                              transform3(point3{7.07f, 7.07f, -250.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(surface_id::e_sensitive,
                              transform3(point3{7.07f, 7.07f, 100.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));

    auto trpz_factory = std::make_shared<trapezoid_factory>();

    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.emplace_back(surface_id::e_sensitive,
                              transform3(point3{7.07f, 7.07f, 600.f}), vol_idx,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_factory->push_back(std::move(trpz_sf_data));

    auto cyl_factory = std::make_shared<cylinder_factory>();

    cyl_sf_data.clear();
    cyl_sf_data.emplace_back(surface_id::e_passive,
                             transform3(point3{0.f, 0.f, 0.f}), vol_idx,
                             std::vector<scalar>{5.f, -1300.f, 1300.f});
    cyl_factory->push_back(std::move(cyl_sf_data));

    // senstivies should end up in the grid, portals in the volume
    gbuilder.add_surfaces(pt_cyl_factory, geo_ctx);
    gbuilder.add_surfaces(rect_factory, geo_ctx);
    gbuilder.add_surfaces(trpz_factory, geo_ctx);
    gbuilder.add_surfaces(cyl_factory, geo_ctx);

    const auto& cyl_axis_z = gbuilder.get().template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.nbins(), 4u);

    gbuilder.build(d);

    //
    // check results
    //
    const auto& vol = d.volumes().back();
    EXPECT_TRUE(d.volumes().size() == 1u);
    EXPECT_EQ(vol.index(), 0u);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);

    // only the portals are referenced through the volume
    typename toy_metadata::object_link_type sf_range{};
    sf_range[0] = {acc_ids::e_default, 0u};
    sf_range[1] = {acc_ids::e_cylinder2_grid, 0u};

    // toy detector makes no distinction between the surface types
    EXPECT_EQ(vol.template accel_link<geo_obj_id::e_portal>(),
              sf_range[geo_obj_id::e_portal]);
    EXPECT_EQ(vol.template accel_link<geo_obj_id::e_sensitive>(),
              sf_range[geo_obj_id::e_sensitive]);
    EXPECT_EQ(vol.template accel_link<geo_obj_id::e_passive>(),
              sf_range[geo_obj_id::e_passive]);

    // Only the portals should be in the detector's surface container now
    EXPECT_EQ(d.surfaces().size(), 7u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1u);

    // check the portals in the detector
    const auto& bf_finder =
        d.accelerator_store()
            .template get<detector_t::accel::id::e_brute_force>()[0];
    for (const auto& sf : bf_finder.all()) {
        EXPECT_TRUE((sf.id() == surface_id::e_portal) or
                    (sf.id() == surface_id::e_passive));
        EXPECT_EQ(sf.volume(), 0u);
    }

    // check the sensitive surfaces in the grid
    const auto& cyl_grid =
        d.accelerator_store()
            .template get<detector_t::accel::id::e_cylinder2_grid>()[0];
    dindex trf_idx{3u};
    for (const auto& sf : cyl_grid.all()) {
        EXPECT_TRUE(sf.is_sensitive());
        EXPECT_EQ(sf.volume(), 0u);
        EXPECT_EQ(sf.transform(), trf_idx++);
    }
}
