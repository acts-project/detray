/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/test/types.hpp"
#include "detray/tools/cuboid_portal_generator.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <memory>

namespace {

using point3 = detray::test::point3;
using vector3 = detray::test::vector3;
using point2 = detray::test::point2;

constexpr detray::scalar tol{1e-7f};

/// Adds a few surfaces to the detector.
template <typename detector_t>
void prefill_detector(detector_t& d,
                      typename detector_t::geometry_context ctx) {
    using scalar_t = typename detector_t::scalar_type;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using surface_t = typename detector_t::surface_type;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    detray::empty_context empty_ctx{};
    vecmem::memory_resource* host_mr = d.resource();
    typename detector_t::transform_container trfs(*host_mr);
    typename detector_t::surface_container_t surfaces{};
    typename detector_t::mask_container masks(*host_mr);
    typename detector_t::material_container materials(*host_mr);

    /// Surface 0
    point3 t0{0.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t0);
    masks.template emplace_back<mask_id::e_rectangle2>(empty_ctx, 0u, -3.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::gold<scalar_t>(), 3.f);
    mask_link_t mask_link{mask_id::e_rectangle2,
                          masks.template size<mask_id::e_rectangle2>() - 1};
    material_link_t material_link{
        material_id::e_slab,
        materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    /// Surface 1
    point3 t1{1.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t1);
    masks.template emplace_back<mask_id::e_annulus2>(empty_ctx, 0u, 1.f, 2.f,
                                                     3.f, 4.f, 5.f, 6.f, 7.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::tungsten<scalar_t>(), 12.f);

    mask_link = {mask_id::e_annulus2,
                 masks.template size<mask_id::e_annulus2>() - 1};
    material_link = {material_id::e_slab,
                     materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    /// Surface 2
    point3 t2{2.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t2);
    masks.template emplace_back<mask_id::e_trapezoid2>(empty_ctx, 0u, 1.f, 2.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_rod>(
        empty_ctx, detray::aluminium<scalar_t>(), 4.f);

    mask_link = {mask_id::e_trapezoid2,
                 masks.template size<mask_id::e_trapezoid2>() - 1};
    material_link = {material_id::e_rod,
                     materials.template size<material_id::e_rod>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    // Add surfaces to lookup, so they can be easily fetched using a barcode
    for (const auto sf : surfaces) {
        d.add_surface_to_lookup(sf);
    }

    // Add the new data
    d.new_volume(detray::volume_id::e_cylinder);
    d.append_portals(std::move(surfaces));
    d.append_transforms(std::move(trfs));
    d.append_masks(std::move(masks));
    d.append_materials(std::move(materials));
}

/// Check volume links for a collection of masks in a given detector
template <typename detector_t, typename detector_t::mask_link::id_type mask_id>
inline void check_mask(const detector_t& d,
                       const std::vector<detray::dindex>& vol_links) {
    for (const auto [idx, mask] :
         detray::views::enumerate(d.mask_store().template get<mask_id>())) {
        EXPECT_EQ(mask.volume_link(), vol_links.at(idx))
            << "mask no. " << idx << ": " << mask.to_string();
    }
}

}  // anonymous namespace

/// This tests the functionality of a detector as a data store manager
GTEST_TEST(detray_core, detector) {

    using namespace detray;

    using detector_t = detector<>;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using finder_id = typename detector_t::sf_finders::id;

    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};

    EXPECT_TRUE(d.volumes().empty());
    EXPECT_TRUE(d.portals().empty());
    EXPECT_TRUE(d.transform_store().empty());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_rectangle2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_rectangle2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_trapezoid2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_annulus2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_cylinder2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_cylinder2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_ring2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_ring2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_straw_wire>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_cell_wire>());
    EXPECT_TRUE(d.material_store().template empty<material_id::e_slab>());
    EXPECT_TRUE(d.material_store().template empty<material_id::e_rod>());
    EXPECT_TRUE(d.surface_store().template empty<finder_id::e_brute_force>());
    EXPECT_TRUE(d.surface_store().template empty<finder_id::e_disc_grid>());
    EXPECT_TRUE(
        d.surface_store().template empty<finder_id::e_cylinder2_grid>());
    EXPECT_TRUE(d.surface_store().template empty<finder_id::e_irr_disc_grid>());
    EXPECT_TRUE(
        d.surface_store().template empty<finder_id::e_irr_cylinder2_grid>());
    EXPECT_TRUE(d.surface_store().template empty<finder_id::e_default>());

    // Add some geometrical data
    prefill_detector(d, geo_ctx);
    // TODO: add B-field check

    EXPECT_EQ(d.volumes().size(), 1u);
    EXPECT_EQ(d.portals().size(), 3u);
    EXPECT_EQ(d.transform_store().size(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_rectangle2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_straw_wire>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cell_wire>(), 0u);
    EXPECT_EQ(d.material_store().template size<material_id::e_slab>(), 2u);
    EXPECT_EQ(d.material_store().template size<material_id::e_rod>(), 1u);
    EXPECT_EQ(d.surface_store().template size<finder_id::e_brute_force>(), 1u);
    EXPECT_EQ(d.surface_store().template size<finder_id::e_disc_grid>(), 0u);
    EXPECT_EQ(d.surface_store().template size<finder_id::e_cylinder2_grid>(),
              0u);
    EXPECT_EQ(d.surface_store().template size<finder_id::e_irr_disc_grid>(),
              0u);
    EXPECT_EQ(
        d.surface_store().template size<finder_id::e_irr_cylinder2_grid>(), 0u);
    EXPECT_EQ(d.surface_store().template size<finder_id::e_default>(), 1u);
}

/// This tests the functionality of a surface factory
GTEST_TEST(detray_tools, surface_factory) {

    using namespace detray;

    using detector_t = detector<>;
    using transform3 = typename detector_t::transform3;
    using mask_id = typename detector_t::masks::id;

    //
    // check portal cylinder
    //
    using portal_cylinder_factory =
        surface_factory<detector_t, typename default_metadata::cylinder_portal,
                        mask_id::e_portal_cylinder2, surface_id::e_portal>;

    auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();

    typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, -1000.f}), 0u,
                             std::vector<scalar>{10.f, -1000.f, 1500.f});
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1000.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1000.f});

    EXPECT_EQ(pt_cyl_factory->size(), 0u);
    EXPECT_TRUE(pt_cyl_factory->components().empty());
    EXPECT_TRUE(pt_cyl_factory->transforms().empty());
    EXPECT_TRUE(pt_cyl_factory->volume_links().empty());

    pt_cyl_factory->push_back(std::move(cyl_sf_data));
    // data should be safely added to the factory by now
    cyl_sf_data.clear();

    EXPECT_EQ(pt_cyl_factory->size(), 2u);
    EXPECT_EQ(pt_cyl_factory->components().size(), 2u);
    EXPECT_EQ(pt_cyl_factory->transforms().size(), 2u);
    EXPECT_EQ(pt_cyl_factory->volume_links().size(), 2u);
    const auto& portal_cyl_comps = pt_cyl_factory->components().front();
    EXPECT_NEAR(portal_cyl_comps[0], 10.f, tol);
    EXPECT_NEAR(portal_cyl_comps[1], -1000.f, tol);
    EXPECT_NEAR(portal_cyl_comps[2], 1500.f, tol);
    const auto& portal_cyl_vol_links = pt_cyl_factory->volume_links();
    EXPECT_EQ(portal_cyl_vol_links[0], 0u);
    EXPECT_EQ(portal_cyl_vol_links[1], 2u);

    //
    // check sensitive cylinder
    //
    using sensitive_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_sensitive>;

    auto sens_cyl_factory = std::make_shared<sensitive_cylinder_factory>();

    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, -50.f}), 1u,
                             std::vector<scalar>{5.f, -900.f, 900.f});
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 50.f}), 1u,
                             std::vector<scalar>{5.f, -900.f, 900.f});

    EXPECT_EQ(sens_cyl_factory->size(), 0u);
    EXPECT_TRUE(sens_cyl_factory->components().empty());
    EXPECT_TRUE(sens_cyl_factory->transforms().empty());
    // The single view is never empty, only uninitialized
    EXPECT_FALSE(sens_cyl_factory->volume_links().empty());

    sens_cyl_factory->push_back(std::move(cyl_sf_data));
    cyl_sf_data.clear();

    EXPECT_EQ(sens_cyl_factory->size(), 2u);
    EXPECT_EQ(sens_cyl_factory->components().size(), 2u);
    EXPECT_EQ(sens_cyl_factory->transforms().size(), 2u);
    EXPECT_EQ(sens_cyl_factory->volume_links().size(), 1u);
    const auto& sens_cyl_comps = sens_cyl_factory->components().front();
    EXPECT_NEAR(sens_cyl_comps[0], 5.f, tol);
    EXPECT_NEAR(sens_cyl_comps[1], -900.f, tol);
    EXPECT_NEAR(sens_cyl_comps[2], 900.f, tol);
    const auto& sens_cyl_vol_links = sens_cyl_factory->volume_links();
    EXPECT_EQ(sens_cyl_vol_links[0], 1u);

    //
    // check passive cylinder
    //
    using passive_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_passive>;

    auto psv_cyl_factory = std::make_shared<passive_cylinder_factory>();

    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, -20.f}), 1u,
                             std::vector<scalar>{4.9f, -900.f, 900.f});
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 20.f}), 1u,
                             std::vector<scalar>{4.9f, -900.f, 900.f});

    EXPECT_EQ(psv_cyl_factory->size(), 0u);
    EXPECT_TRUE(psv_cyl_factory->components().empty());
    EXPECT_TRUE(psv_cyl_factory->transforms().empty());
    // The single view is never empty, only uninitialized
    EXPECT_FALSE(psv_cyl_factory->volume_links().empty());

    psv_cyl_factory->push_back(std::move(cyl_sf_data));
    cyl_sf_data.clear();

    EXPECT_EQ(psv_cyl_factory->size(), 2u);
    EXPECT_EQ(psv_cyl_factory->components().size(), 2u);
    EXPECT_EQ(psv_cyl_factory->transforms().size(), 2u);
    EXPECT_EQ(psv_cyl_factory->volume_links().size(), 1u);
    const auto& psv_cyl_comps = psv_cyl_factory->components().front();
    EXPECT_NEAR(psv_cyl_comps[0], 4.9f, tol);
    EXPECT_NEAR(psv_cyl_comps[1], -900.f, tol);
    EXPECT_NEAR(psv_cyl_comps[2], 900.f, tol);
    const auto& psv_cyl_vol_links = psv_cyl_factory->volume_links();
    EXPECT_EQ(psv_cyl_vol_links[0], 1u);

    //
    // check the other mask types for this detector
    //

    // annulus
    using annulus_factory =
        surface_factory<detector_t, annulus2D<>, mask_id::e_annulus2,
                        surface_id::e_sensitive>;

    auto ann_factory = std::make_shared<annulus_factory>();

    typename annulus_factory::sf_data_collection ann_sf_data;
    ann_sf_data.emplace_back(
        transform3(point3{0.f, 0.f, 0.f}), 1u,
        std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
    ann_factory->push_back(std::move(ann_sf_data));
    ann_sf_data.clear();

    const auto& ann_comps = ann_factory->components().front();
    EXPECT_NEAR(ann_comps[0], 300.f, tol);
    EXPECT_NEAR(ann_comps[1], 350.f, tol);
    EXPECT_NEAR(ann_comps[2], -0.1f, tol);
    EXPECT_NEAR(ann_comps[3], 0.1f, tol);
    EXPECT_NEAR(ann_comps[4], 0.5f, tol);
    EXPECT_NEAR(ann_comps[5], 0.6f, tol);
    EXPECT_NEAR(ann_comps[6], 1.4f, tol);

    // rectangles
    using rectangle_factory =
        surface_factory<detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;

    auto rect_factory = std::make_shared<rectangle_factory>();

    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 1u,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));
    rect_sf_data.clear();

    const auto& rectgl_comps = rect_factory->components().front();
    EXPECT_NEAR(rectgl_comps[0], 10.f, tol);
    EXPECT_NEAR(rectgl_comps[1], 8.f, tol);

    // ring
    using ring_factory = surface_factory<detector_t, ring2D<>, mask_id::e_ring2,
                                         surface_id::e_passive>;

    auto rng_factory = std::make_shared<ring_factory>();

    typename ring_factory::sf_data_collection ring_sf_data;
    ring_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 1u,
                              std::vector<scalar>{0.f, 5.f});
    rng_factory->push_back(std::move(ring_sf_data));
    ring_sf_data.clear();

    const auto& ring_comps = rng_factory->components().front();
    EXPECT_NEAR(ring_comps[0], 0.f, tol);
    EXPECT_NEAR(ring_comps[1], 5.f, tol);

    // trapezoid
    using trapezoid_factory =
        surface_factory<detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;

    auto trpz_factory = std::make_shared<trapezoid_factory>();

    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 1u,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_factory->push_back(std::move(trpz_sf_data));
    trpz_sf_data.clear();

    const auto& trpz_comps = trpz_factory->components().front();
    EXPECT_NEAR(trpz_comps[0], 1.f, tol);
    EXPECT_NEAR(trpz_comps[1], 3.f, tol);
    EXPECT_NEAR(trpz_comps[2], 2.f, tol);
    EXPECT_NEAR(trpz_comps[3], 0.25f, tol);
}

/// This tests the initialization of a detector volume using a volume builder
GTEST_TEST(detray_tools, volume_builder) {

    using namespace detray;

    vecmem::host_memory_resource host_mr;

    using detector_t = detector<>;

    detector_t d(host_mr);

    EXPECT_TRUE(d.volumes().size() == 0u);

    volume_builder<detector_t> vbuilder{volume_id::e_cylinder};
    vbuilder.build(d);

    const auto& vol = d.volumes().back();

    EXPECT_TRUE(d.volumes().size() == 1u);
    EXPECT_EQ(vol.index(), 0u);
    EXPECT_EQ(vol.index(), vbuilder.vol_index());
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);
}

/// Integration test to build a cylinder volume with contained surfaces
GTEST_TEST(detray_tools, detector_volume_construction) {

    using namespace detray;

    using detector_t = detector<>;
    using transform3 = typename detector_t::transform3;
    using geo_obj_id = typename detector_t::geo_obj_ids;
    using mask_id = typename detector_t::masks::id;
    using sf_finder_id = typename detector_t::sf_finders::id;

    // portal factories
    using portal_cylinder_factory =
        surface_factory<detector_t, typename default_metadata::cylinder_portal,
                        mask_id::e_portal_cylinder2, surface_id::e_portal>;
    using portal_disc_factory =
        surface_factory<detector_t, ring2D<>, mask_id::e_portal_ring2,
                        surface_id::e_portal>;

    // sensitive/passive surface factories
    using annulus_factory =
        surface_factory<detector_t, annulus2D<>, mask_id::e_annulus2,
                        surface_id::e_sensitive>;
    using cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_passive>;
    using rectangle_factory =
        surface_factory<detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;
    using ring_factory = surface_factory<detector_t, ring2D<>, mask_id::e_ring2,
                                         surface_id::e_passive>;
    using trapezoid_factory =
        surface_factory<detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;

    // detector
    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};
    // ensure there is a data offset that needs to be handled correctly
    prefill_detector(d, geo_ctx);
    const dindex first_trf{d.transform_store().size()};
    const auto vol_idx{
        static_cast<typename detector_t::surface_type::navigation_link>(
            d.volumes().size())};

    // initial checks
    EXPECT_EQ(d.volumes().size(), 1u);
    EXPECT_EQ(d.portals().size(), 3u);

    // volume builder
    volume_builder<detector_t> vbuilder{volume_id::e_cylinder};
    typename detector_t::point3 t{0.f, 0.f, 20.f};
    vbuilder.add_volume_placement(t);

    //
    // Fill the surface factories with data
    //

    // portal surfaces
    auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();
    typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
    // Creates two cylinders at radius 0mm and 10mm with an extent in z
    // of -5mm to 5mm and linking to volumes 0 and 2, respectively.
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 0u,
                             std::vector<scalar>{10.f, -1500.f, 1500.f});
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1500.f});
    pt_cyl_factory->push_back(std::move(cyl_sf_data));

    auto pt_disc_factory = std::make_shared<portal_disc_factory>();
    typename portal_disc_factory::sf_data_collection ring_sf_data;
    // Creates two discs with a radius of 10mm, linking to volumes 3 and 4
    ring_sf_data.emplace_back(transform3(point3{0.f, 0.f, -1500.f}), 3u,
                              std::vector<scalar>{0.f, 10.f});
    ring_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1500.f}), 4u,
                              std::vector<scalar>{0.f, 10.f});
    pt_disc_factory->push_back(std::move(ring_sf_data));

    // sensitive surfaces
    auto ann_factory = std::make_shared<annulus_factory>();
    typename annulus_factory::sf_data_collection ann_sf_data;
    ann_sf_data.emplace_back(
        transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
        std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
    ann_sf_data.emplace_back(
        transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
        std::vector<scalar>{350.f, 400.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
    ann_factory->push_back(std::move(ann_sf_data));

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

    // passive surfaces
    auto cyl_factory = std::make_shared<cylinder_factory>();
    cyl_sf_data.clear();
    cyl_sf_data.emplace_back(transform3(point3{0.f, 0.f, 0.f}), vol_idx,
                             std::vector<scalar>{5.f, -1300.f, 1300.f});
    cyl_factory->push_back(std::move(cyl_sf_data));

    auto rng_factory = std::make_shared<ring_factory>();
    ring_sf_data.clear();
    ring_sf_data.emplace_back(transform3(point3{0.f, 0.f, -1300.f}), vol_idx,
                              std::vector<scalar>{0.f, 5.f});
    ring_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1300.f}), vol_idx,
                              std::vector<scalar>{0.f, 5.f});
    rng_factory->push_back(std::move(ring_sf_data));

    //
    // Fill everything into volume
    //
    vbuilder.add_portals(pt_cyl_factory, geo_ctx);
    vbuilder.add_portals(pt_disc_factory, geo_ctx);

    vbuilder.add_sensitives(ann_factory, geo_ctx);
    vbuilder.add_sensitives(rect_factory, geo_ctx);
    vbuilder.add_sensitives(trpz_factory, geo_ctx);

    vbuilder.add_passives(cyl_factory, geo_ctx);
    vbuilder.add_passives(rng_factory, geo_ctx);

    // try adding something extra later...
    rect_factory->clear();
    rect_sf_data.clear();
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, 10.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, 20.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));
    rect_sf_data.clear();
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, 30.f}), vol_idx,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));

    vbuilder.add_sensitives(rect_factory, geo_ctx);

    //
    // Adds all surfaces to the detector
    //
    vbuilder.build(d);

    //
    // check results
    //
    const auto& vol = d.volumes().back();

    EXPECT_EQ(d.volumes().size(), 2u);
    EXPECT_EQ(vol.index(), 1u);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);

    // Check the volume placement
    typename detector_t::transform3 trf{t};
    EXPECT_TRUE(d.transform_store()[first_trf] == trf);

    // Check the acceleration data structure link
    dtyped_index<sf_finder_id, dindex> acc_link{sf_finder_id::e_default, 1u};
    ASSERT_TRUE(vol.full_link().size() == geo_obj_id::e_size);
    EXPECT_EQ(vol.link<geo_obj_id::e_portal>(), acc_link);
    EXPECT_EQ(vol.link<geo_obj_id::e_passive>(), acc_link);
    // Not set by the vanilla volume builder
    EXPECT_TRUE(is_invalid_value(vol.link<geo_obj_id::e_sensitive>()));

    EXPECT_EQ(d.portals().size(), 19u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 2u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 4u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 7u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 4u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 2u);

    // check surface type and volume link
    std::vector<surface_id> sf_ids{};
    sf_ids.reserve(d.n_surfaces());
    sf_ids.insert(sf_ids.end(), 3u, surface_id::e_sensitive);
    sf_ids.insert(sf_ids.end(), 4u, surface_id::e_portal);
    sf_ids.insert(sf_ids.end(), 6u, surface_id::e_sensitive);
    sf_ids.insert(sf_ids.end(), 3u, surface_id::e_passive);
    sf_ids.insert(sf_ids.end(), 3u, surface_id::e_sensitive);

    std::vector<dindex> volume_links{};
    volume_links.reserve(d.n_surfaces());
    volume_links.insert(volume_links.end(), 3u, 0u);
    volume_links.insert(volume_links.end(), 16u, 1u);

    // Check surface id and volume links
    for (const auto [idx, sf_id] : detray::views::enumerate(sf_ids)) {
        geometry::barcode bcd{};
        bcd.set_index(idx);
        const auto& sf = d.surface(bcd);
        EXPECT_EQ(sf.id(), sf_id) << "error at index: " << idx;
        EXPECT_EQ(sf.volume(), volume_links.at(idx))
            << "error at index: " << idx;
    }

    // check that the transform indices are continuous for the newly added
    // surfaces. The first new transform belongs to the volume itself
    for (std::size_t idx :
         detray::views::iota(dindex_range{3, d.n_surfaces()})) {
        geometry::barcode bcd{};
        bcd.set_index(idx);
        EXPECT_EQ(d.surface(bcd).transform(), idx + 1)
            << "error at index: " << idx;
    }

    // check surface mask links
    std::vector<typename detector_t::mask_link> mask_links{
        {mask_id::e_rectangle2, 0u},
        {mask_id::e_annulus2, 0u},
        {mask_id::e_trapezoid2, 0u},
        {mask_id::e_portal_cylinder2, 0u},
        {mask_id::e_portal_cylinder2, 1u},
        {mask_id::e_portal_ring2, 0u},
        {mask_id::e_portal_ring2, 1u},
        {mask_id::e_annulus2, 1u},
        {mask_id::e_annulus2, 2u},
        {mask_id::e_rectangle2, 1u},
        {mask_id::e_rectangle2, 2u},
        {mask_id::e_rectangle2, 3u},
        {mask_id::e_trapezoid2, 1u},
        {mask_id::e_cylinder2, 0u},
        {mask_id::e_ring2, 2u},
        {mask_id::e_ring2, 3u},
        {mask_id::e_rectangle2, 4u},
        {mask_id::e_rectangle2, 5u},
        {mask_id::e_rectangle2, 6u}};
    for (const auto [idx, m_link] : detray::views::enumerate(mask_links)) {
        geometry::barcode bcd{};
        bcd.set_index(idx);
        EXPECT_EQ(d.surface(bcd).mask(), m_link) << "error at index: " << idx;
    }

    // check mask volume links
    volume_links.clear();
    volume_links = {0u, 2u};
    check_mask<detector_t, mask_id::e_portal_cylinder2>(d, volume_links);

    volume_links.clear();
    volume_links = {1u};
    check_mask<detector_t, mask_id::e_cylinder2>(d, volume_links);

    volume_links.clear();
    volume_links = {3u, 4u, 1u, 1u};
    check_mask<detector_t, mask_id::e_portal_ring2>(d, volume_links);
    check_mask<detector_t, mask_id::e_ring2>(d, volume_links);

    volume_links.clear();
    volume_links = {0u, 1u, 1u};
    check_mask<detector_t, mask_id::e_annulus2>(d, volume_links);

    volume_links.clear();
    volume_links.reserve(7u);
    volume_links.push_back(0u);
    volume_links.insert(volume_links.end(), 6u, 1u);
    check_mask<detector_t, mask_id::e_rectangle2>(d, volume_links);

    volume_links.clear();
    volume_links = {0u, 1u};
    check_mask<detector_t, mask_id::e_trapezoid2>(d, volume_links);
}

/// Integration test to build a cylinder volume with contained surfaces
GTEST_TEST(detray_tools, detector_builder) {
    using namespace detray;

    using detector_t = detector<>;
    using transform3 = typename detector_t::transform3;
    using mask_id = typename detector_t::masks::id;

    // Surface factories
    using trapezoid_factory =
        surface_factory<detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;

    // detector builder
    auto geo_ctx = typename detector_t::geometry_context{};

    detector_builder<default_metadata> det_builder{};

    //
    // first volume builder
    //
    auto vbuilder = det_builder.new_volume(volume_id::e_cylinder);

    typename detector_t::point3 t{0.f, 0.f, 0.f};
    vbuilder->add_volume_placement(t);

    // initial checks
    EXPECT_EQ(vbuilder->vol_index(), 0u);

    const auto vol_idx{
        static_cast<typename detector_t::surface_type::navigation_link>(
            vbuilder->vol_index())};

    auto trpz_factory = std::make_shared<trapezoid_factory>();
    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1100.f}), vol_idx,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_sf_data.emplace_back(transform3(point3{0.f, 0.f, 1200.f}), vol_idx,
                              std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
    trpz_factory->push_back(std::move(trpz_sf_data));

    vbuilder->add_sensitives(trpz_factory, geo_ctx);

    //
    // second volume builder
    //
    auto vbuilder2 = det_builder.new_volume(volume_id::e_cuboid);

    // volume builder
    t = typename detector_t::point3{0.f, 0.f, 20.f};
    vbuilder2->add_volume_placement(t);

    // Add a portal box around the cuboid volume with a min distance of 'env'
    constexpr auto env{0.1f * unit<detray::scalar>::mm};
    auto portal_generator =
        std::make_shared<cuboid_portal_generator<detector_t>>(env);

    vbuilder2->add_portals(portal_generator);

    // initial checks
    EXPECT_EQ(vbuilder2->vol_index(), 1u);

    //
    // build the detector
    //
    vecmem::host_memory_resource host_mr;
    const detector_t d = det_builder.build(host_mr);
    const auto& vol0 = d.volume_by_index(0u);
    const auto& vol1 = d.volume_by_index(1u);

    // check the results
    EXPECT_EQ(d.volumes().size(), 2u);
    EXPECT_EQ(vol0.id(), volume_id::e_cylinder);
    EXPECT_EQ(vol0.index(), 0u);
    EXPECT_EQ(vol1.id(), volume_id::e_cuboid);
    EXPECT_EQ(vol1.index(), 1u);

    // Check the volume placements for both volumes
    typename detector_t::transform3 identity{};
    EXPECT_TRUE(vol0.transform() == identity);
    EXPECT_TRUE(d.transform_store()[0u] == identity);
    typename detector_t::transform3 trf{t};
    EXPECT_TRUE(vol1.transform() == trf);
    EXPECT_TRUE(d.transform_store()[4u] == trf);

    // Check the acceleration data structure link (indirectly)
    EXPECT_EQ(vol0.n_max_candidates(), 3u);
    EXPECT_EQ(vol1.n_max_candidates(), 6u);

    EXPECT_EQ(d.surface_lookup().size(), 9u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 3u);
}
