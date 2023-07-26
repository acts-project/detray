/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"
#include "detray/test/types.hpp"
#include "detray/tools/grid_builder.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <limits>

using namespace detray;
using namespace detray::n_axis;

namespace {

using point3 = test::point3;
using vector3 = test::vector3;

test::transform3 Identity{};

using test_detector_t = detector<toy_metadata<>>;

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_tools, grid_factory) {

    // Data-owning grid collection
    vecmem::host_memory_resource host_mr;
    auto gr_factory =
        grid_factory<dindex, simple_serializer, regular_attacher<3>>{host_mr};

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
    auto loc_p = ann_gr.global_to_local(Identity, p, d);
    ann_gr.populate(loc_p, 3u);

    EXPECT_EQ(ann_gr.search(loc_p)[0], 3u);

    // gr_builder.to_string(ann_gr);

    // Build from parameters
    const std::vector<scalar> bin_edges_z{-10.f, -8.f, -6.5f, -1.f,
                                          4.f,   5.f,  6.f,   9.f};
    const std::vector<scalar> bin_edges_phi{};

    auto cyl_gr = gr_factory.template new_grid<cylinder2D<>>(
        {bin_edges_z.front(), bin_edges_z.back(), 0.f,
         2.f * constant<scalar>::pi},
        {bin_edges_z.size() - 1, 10u}, {bin_edges_phi, bin_edges_z},
        std::tuple<circular<label::e_rphi>, closed<label::e_cyl_z>>{},
        std::tuple<regular<>, irregular<>>{});

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
    loc_p = cyl_gr.global_to_local(Identity, p, d);
    cyl_gr.populate(loc_p, 33u);

    EXPECT_EQ(cyl_gr.search(loc_p)[0], 33u);

    // Build the same cylinder grid from a mask
    const scalar r{5.f};
    const scalar n_half_z{-10.f};
    const scalar p_half_z{9.f};
    mask<cylinder2D<>> cyl2{0u, r, n_half_z, p_half_z};

    auto cyl_gr2 = gr_factory.template new_grid<circular<label::e_rphi>,
                                                closed<label::e_cyl_z>,
                                                regular<>, irregular<>>(
        cyl2, {bin_edges_z.size() - 1, 10u}, {bin_edges_phi, bin_edges_z});

    // Get grid collection with the correct allocator
    auto grid_coll = gr_factory.new_collection<decltype(cyl_gr)>();
    EXPECT_TRUE(grid_coll.size() == 0u);
    grid_coll.push_back(cyl_gr);
    grid_coll.push_back(cyl_gr2);
    EXPECT_TRUE(grid_coll.size() == 2u);

    EXPECT_EQ(grid_coll[0].search(loc_p)[0], 33u);
    // gr_factory.to_string(grid_coll[0]);
}

/// Unittest: Test the grid builder
GTEST_TEST(detray_tools, grid_builder) {

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             test_detector_t::surface_type, simple_serializer,
             regular_attacher<9>>;

    auto gbuilder = grid_builder<test_detector_t, cyl_grid_t,
                                 detray::detail::bin_associator>{};

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

    using transform3 = typename test_detector_t::transform3;
    using geo_obj_id = typename test_detector_t::geo_obj_ids;
    using acc_ids = typename test_detector_t::sf_finders::id;
    using mask_id = typename test_detector_t::masks::id;

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             test_detector_t::surface_type, simple_serializer,
             regular_attacher<9>>;

    using portal_cylinder_factory_t =
        surface_factory<test_detector_t, cylinder2D<>,
                        mask_id::e_portal_cylinder2, surface_id::e_portal>;
    using rectangle_factory =
        surface_factory<test_detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;
    using trapezoid_factory =
        surface_factory<test_detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;
    using cylinder_factory =
        surface_factory<test_detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_passive>;

    vecmem::host_memory_resource host_mr;
    test_detector_t d(host_mr);
    auto geo_ctx = typename test_detector_t::geometry_context{};
    const auto vol_idx{
        static_cast<typename test_detector_t::surface_type::navigation_link>(
            d.volumes().size())};

    auto vbuilder = std::make_unique<volume_builder<test_detector_t>>(
        volume_id::e_cylinder);
    auto gbuilder =
        grid_builder<test_detector_t, cyl_grid_t>{std::move(vbuilder)};
    // passive surfaces are added to the grid
    // gbuilder.set_add_passives();

    // The cylinder portals are at the end of the surface range by construction
    const auto cyl_mask = mask<cylinder2D<>>{0u, 10.f, -500.f, 500.f};
    std::size_t n_phi_bins{5u}, n_z_bins{4u};

    // Build empty grid
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});

    EXPECT_TRUE(d.volumes().size() == 0);

    // Add some portals first
    auto pt_cyl_factory = std::make_shared<portal_cylinder_factory_t>();

    typename portal_cylinder_factory_t::sf_data_collection cyl_sf_data;
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 0u,
                             std::vector<scalar>{10.f, -1500.f, 1500.f});
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1500.f});
    pt_cyl_factory->push_back(std::move(cyl_sf_data));

    // Then some passive and sensitive surfaces
    auto rect_factory = std::make_shared<rectangle_factory>();

    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -10.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -20.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -30.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));

    auto trpz_factory = std::make_shared<trapezoid_factory>();

    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_factory->push_back(std::move(trpz_sf_data));

    auto cyl_factory = std::make_shared<cylinder_factory>();

    cyl_sf_data.clear();
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), vol_idx,
                             std::vector<scalar>{5.f, -1300.f, 1300.f});
    cyl_factory->push_back(std::move(cyl_sf_data));

    // senstivies should end up in the grid, portals in the volume
    gbuilder.add_portals(pt_cyl_factory, geo_ctx);
    gbuilder.add_sensitives(rect_factory, geo_ctx);
    gbuilder.add_sensitives(trpz_factory, geo_ctx);
    gbuilder.add_passives(cyl_factory, geo_ctx);

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
    typename toy_metadata<>::object_link_type sf_range{};
    sf_range[0] = {acc_ids::e_default, 0u};
    sf_range[1] = {acc_ids::e_cylinder2_grid, 0u};
    // toy detector makes no distinction between the surface types
    EXPECT_EQ(vol.template link<geo_obj_id::e_portal>(),
              sf_range[geo_obj_id::e_portal]);
    EXPECT_EQ(vol.template link<geo_obj_id::e_sensitive>(),
              sf_range[geo_obj_id::e_sensitive]);
    EXPECT_EQ(vol.template link<geo_obj_id::e_passive>(),
              sf_range[geo_obj_id::e_passive]);

    // Only the portals should be in the detector's surface container now
    EXPECT_EQ(d.surface_lookup().size(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1u);

    // check the portals in the detector
    for (const auto& sf : d.surface_lookup()) {
        EXPECT_TRUE((sf.id() == surface_id::e_portal) or
                    (sf.id() == surface_id::e_passive));
    }

    // check the sensitive surfaces in the grid
    const auto& cyl_grid =
        d.surface_store()
            .template get<
                test_detector_t::sf_finders::id::e_cylinder2_grid>()[0];
    dindex trf_idx{4u};
    for (const auto& sf : cyl_grid.all()) {
        EXPECT_TRUE(sf.is_sensitive());
        EXPECT_EQ(sf.volume(), 0u);
        EXPECT_EQ(sf.transform(), trf_idx++);
    }
}
