/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <memory>

namespace {

/// @note __plugin has to be defined with a preprocessor command
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

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
    typename detector_t::surface_container surfaces{};
    typename detector_t::mask_container masks(*host_mr);
    typename detector_t::material_container materials(*host_mr);

    /// Surface 0
    point3 t0{0.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t0);
    masks.template emplace_back<mask_id::e_rectangle2>(empty_ctx, 0UL, -3.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::gold<scalar_t>(), 3.f);
    mask_link_t mask_link{mask_id::e_rectangle2,
                          masks.template size<mask_id::e_rectangle2>() - 1};
    material_link_t material_link{
        material_id::e_slab,
        materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0UL,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);
    /// Surface 1
    point3 t1{1.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t1);
    masks.template emplace_back<mask_id::e_annulus2>(empty_ctx, 0UL, 1.f, 2.f,
                                                     3.f, 4.f, 5.f, 6.f, 7.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::tungsten<scalar_t>(), 12.f);

    mask_link = {mask_id::e_annulus2,
                 masks.template size<mask_id::e_annulus2>() - 1};
    material_link = {material_id::e_slab,
                     materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0UL,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);

    /// Surface 2
    point3 t2{2.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t2);
    masks.template emplace_back<mask_id::e_trapezoid2>(empty_ctx, 0UL, 1.f, 2.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_rod>(
        empty_ctx, detray::aluminium<scalar_t>(), 4.f);

    mask_link = {mask_id::e_trapezoid2,
                 masks.template size<mask_id::e_trapezoid2>() - 1};
    material_link = {material_id::e_rod,
                     materials.template size<material_id::e_rod>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0UL,
                          detray::dindex_invalid,
                          detray::surface_id::e_sensitive);

    // Add the new data
    d.new_volume(detray::volume_id::e_cylinder,
                 {0.f, 10.f, -5.f, 5.f, -M_PI, M_PI});
    d.append_surfaces(std::move(surfaces));
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
TEST(detector, detector) {

    using namespace detray;

    using detector_t =
        detector<detector_registry::default_detector, covfie::field>;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using finder_id = typename detector_t::sf_finders::id;

    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};

    EXPECT_TRUE(d.volumes().empty());
    EXPECT_TRUE(d.surfaces().empty());
    EXPECT_TRUE(d.transform_store().empty());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_rectangle2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_trapezoid2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_annulus2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_cylinder2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_cylinder2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_ring2>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_ring2>());
    EXPECT_TRUE(d.material_store().template empty<material_id::e_slab>());
    EXPECT_TRUE(d.material_store().template empty<material_id::e_rod>());
    EXPECT_TRUE(d.sf_finder_store().template empty<finder_id::e_brute_force>());
    EXPECT_TRUE(d.sf_finder_store().template empty<finder_id::e_disc_grid>());
    EXPECT_TRUE(
        d.sf_finder_store().template empty<finder_id::e_cylinder_grid>());
    EXPECT_TRUE(d.sf_finder_store().template empty<finder_id::e_default>());

    // Add some geometrical data
    prefill_detector(d, geo_ctx);
    // TODO: add B-field check

    EXPECT_EQ(d.volumes().size(), 1UL);
    EXPECT_EQ(d.surfaces().size(), 3UL);
    EXPECT_EQ(d.transform_store().size(), 3UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 1UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 1UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 0UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 0UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 0UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0UL);
    EXPECT_EQ(d.material_store().template size<material_id::e_slab>(), 2UL);
    EXPECT_EQ(d.material_store().template size<material_id::e_rod>(), 1UL);
    EXPECT_EQ(d.sf_finder_store().template size<finder_id::e_brute_force>(),
              0UL);
    EXPECT_EQ(d.sf_finder_store().template size<finder_id::e_disc_grid>(), 0UL);
    EXPECT_EQ(d.sf_finder_store().template size<finder_id::e_cylinder_grid>(),
              0UL);
    EXPECT_EQ(d.sf_finder_store().template size<finder_id::e_default>(), 0UL);
}

/// This tests the functionality of a surface factory
TEST(detector, surface_factory) {

    using namespace detray;

    using detector_t =
        detector<detector_registry::default_detector, covfie::field>;
    using transform3 = typename detector_t::transform3;
    using mask_id = typename detector_t::masks::id;

    //
    // check portal cylinder
    //
    using portal_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_portal_cylinder2,
                        surface_id::e_portal>;

    auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();

    typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, -1000.f}), 0UL,
            std::vector<scalar>{10.f, -1000.f, 1500.f}));
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 1000.f}), 2UL,
            std::vector<scalar>{20.f, -1500.f, 1000.f}));

    EXPECT_EQ(pt_cyl_factory->size(), 0UL);
    EXPECT_TRUE(pt_cyl_factory->components().empty());
    EXPECT_TRUE(pt_cyl_factory->transforms().empty());
    EXPECT_TRUE(pt_cyl_factory->volume_links().empty());

    pt_cyl_factory->add_components(std::move(cyl_sf_data));
    // data should be safely added to the factory by now
    cyl_sf_data.clear();

    EXPECT_EQ(pt_cyl_factory->size(), 2UL);
    EXPECT_EQ(pt_cyl_factory->components().size(), 2UL);
    EXPECT_EQ(pt_cyl_factory->transforms().size(), 2UL);
    EXPECT_EQ(pt_cyl_factory->volume_links().size(), 2UL);
    const auto& portal_cyl_comps = pt_cyl_factory->components().front();
    EXPECT_FLOAT_EQ(portal_cyl_comps[0], scalar{10.f});
    EXPECT_FLOAT_EQ(portal_cyl_comps[1], scalar{-1000.f});
    EXPECT_FLOAT_EQ(portal_cyl_comps[2], scalar{1500.f});
    const auto& portal_cyl_vol_links = pt_cyl_factory->volume_links();
    EXPECT_EQ(portal_cyl_vol_links[0], 0UL);
    EXPECT_EQ(portal_cyl_vol_links[1], 2UL);

    //
    // check sensitive cylinder
    //
    using sensitive_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_sensitive>;

    auto sens_cyl_factory = std::make_shared<sensitive_cylinder_factory>();

    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, -50.f}), 1UL,
            std::vector<scalar>{5.f, -900.f, 900.f}));
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 50.f}), 1UL,
            std::vector<scalar>{5.f, -900.f, 900.f}));

    EXPECT_EQ(sens_cyl_factory->size(), 0UL);
    EXPECT_TRUE(sens_cyl_factory->components().empty());
    EXPECT_TRUE(sens_cyl_factory->transforms().empty());
    // The single view is never empty, only uninitialized
    EXPECT_FALSE(sens_cyl_factory->volume_links().empty());

    sens_cyl_factory->add_components(std::move(cyl_sf_data));
    cyl_sf_data.clear();

    EXPECT_EQ(sens_cyl_factory->size(), 2UL);
    EXPECT_EQ(sens_cyl_factory->components().size(), 2UL);
    EXPECT_EQ(sens_cyl_factory->transforms().size(), 2UL);
    EXPECT_EQ(sens_cyl_factory->volume_links().size(), 1UL);
    const auto& sens_cyl_comps = sens_cyl_factory->components().front();
    EXPECT_FLOAT_EQ(sens_cyl_comps[0], scalar{5.f});
    EXPECT_FLOAT_EQ(sens_cyl_comps[1], scalar{-900.f});
    EXPECT_FLOAT_EQ(sens_cyl_comps[2], scalar{900.f});
    const auto& sens_cyl_vol_links = sens_cyl_factory->volume_links();
    EXPECT_EQ(sens_cyl_vol_links[0], 1UL);

    //
    // check passive cylinder
    //
    using passive_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_passive>;

    auto psv_cyl_factory = std::make_shared<passive_cylinder_factory>();

    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, -20.f}), 1UL,
            std::vector<scalar>{4.9f, -900.f, 900.f}));
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 20.f}), 1UL,
            std::vector<scalar>{4.9f, -900.f, 900.f}));

    EXPECT_EQ(psv_cyl_factory->size(), 0UL);
    EXPECT_TRUE(psv_cyl_factory->components().empty());
    EXPECT_TRUE(psv_cyl_factory->transforms().empty());
    // The single view is never empty, only uninitialized
    EXPECT_FALSE(psv_cyl_factory->volume_links().empty());

    psv_cyl_factory->add_components(std::move(cyl_sf_data));
    cyl_sf_data.clear();

    EXPECT_EQ(psv_cyl_factory->size(), 2UL);
    EXPECT_EQ(psv_cyl_factory->components().size(), 2UL);
    EXPECT_EQ(psv_cyl_factory->transforms().size(), 2UL);
    EXPECT_EQ(psv_cyl_factory->volume_links().size(), 1UL);
    const auto& psv_cyl_comps = psv_cyl_factory->components().front();
    EXPECT_FLOAT_EQ(psv_cyl_comps[0], scalar{4.9f});
    EXPECT_FLOAT_EQ(psv_cyl_comps[1], scalar{-900.f});
    EXPECT_FLOAT_EQ(psv_cyl_comps[2], scalar{900.f});
    const auto& psv_cyl_vol_links = psv_cyl_factory->volume_links();
    EXPECT_EQ(psv_cyl_vol_links[0], 1UL);

    //
    // check the other mask types for this detector
    //

    // annulus
    using annulus_factory =
        surface_factory<detector_t, annulus2D<>, mask_id::e_annulus2,
                        surface_id::e_sensitive>;

    auto ann_factory = std::make_shared<annulus_factory>();

    typename annulus_factory::sf_data_collection ann_sf_data;
    ann_sf_data.push_back(
        std::make_unique<surface_data<detector_t, annulus2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), 1UL,
            std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f}));
    ann_factory->add_components(std::move(ann_sf_data));
    ann_sf_data.clear();

    const auto& ann_comps = ann_factory->components().front();
    EXPECT_FLOAT_EQ(ann_comps[0], scalar{300.f});
    EXPECT_FLOAT_EQ(ann_comps[1], scalar{350.f});
    EXPECT_FLOAT_EQ(ann_comps[2], scalar{-0.1f});
    EXPECT_FLOAT_EQ(ann_comps[3], scalar{0.1f});
    EXPECT_FLOAT_EQ(ann_comps[4], scalar{0.5f});
    EXPECT_FLOAT_EQ(ann_comps[5], scalar{0.6f});
    EXPECT_FLOAT_EQ(ann_comps[6], scalar{1.4f});

    // rectangles
    using rectangle_factory =
        surface_factory<detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;

    auto rect_factory = std::make_shared<rectangle_factory>();

    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), 1UL,
            std::vector<scalar>{10.f, 8.f}));
    rect_factory->add_components(std::move(rect_sf_data));
    rect_sf_data.clear();

    const auto& rectgl_comps = rect_factory->components().front();
    EXPECT_FLOAT_EQ(rectgl_comps[0], scalar{10.f});
    EXPECT_FLOAT_EQ(rectgl_comps[1], scalar{8.f});

    // ring
    using ring_factory = surface_factory<detector_t, ring2D<>, mask_id::e_ring2,
                                         surface_id::e_passive>;

    auto rng_factory = std::make_shared<ring_factory>();

    typename ring_factory::sf_data_collection ring_sf_data;
    ring_sf_data.push_back(std::make_unique<surface_data<detector_t, ring2D<>>>(
        transform3(point3{0.f, 0.f, 0.f}), 1UL, std::vector<scalar>{0.f, 5.f}));
    rng_factory->add_components(std::move(ring_sf_data));
    ring_sf_data.clear();

    const auto& ring_comps = rng_factory->components().front();
    EXPECT_FLOAT_EQ(ring_comps[0], scalar{0.f});
    EXPECT_FLOAT_EQ(ring_comps[1], scalar{5.f});

    // trapezoid
    using trapezoid_factory =
        surface_factory<detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;

    auto trpz_factory = std::make_shared<trapezoid_factory>();

    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.push_back(
        std::make_unique<surface_data<detector_t, trapezoid2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), 1UL,
            std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}));
    trpz_factory->add_components(std::move(trpz_sf_data));
    trpz_sf_data.clear();

    const auto& trpz_comps = trpz_factory->components().front();
    EXPECT_FLOAT_EQ(trpz_comps[0], scalar{1.f});
    EXPECT_FLOAT_EQ(trpz_comps[1], scalar{3.f});
    EXPECT_FLOAT_EQ(trpz_comps[2], scalar{2.f});
    EXPECT_FLOAT_EQ(trpz_comps[3], scalar{0.25f});
}

/// This tests the initialization of a detector volume using a volume builder
TEST(detector, volume_builder) {

    using namespace detray;

    vecmem::host_memory_resource host_mr;

    using detector_t =
        detector<detector_registry::default_detector, covfie::field>;

    detector_t d(host_mr);

    EXPECT_TRUE(d.volumes().size() == 0UL);

    volume_builder<detector_t> vbuilder{};
    const std::array<scalar, 6UL> bounds{0.f, 10.f, -5.f, 5.f, -M_PI, M_PI};
    vbuilder.init_vol(d, volume_id::e_cylinder, bounds);
    const auto& vol = d.volumes().back();

    EXPECT_TRUE(d.volumes().size() == 1UL);
    EXPECT_TRUE(vol.empty());
    EXPECT_EQ(vol.index(), 0UL);
    EXPECT_EQ(vol.index(), vbuilder.get_vol_index());
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);
    for (const auto [idx, b] : detray::views::enumerate(vol.bounds())) {
        EXPECT_FLOAT_EQ(b, bounds[idx]) << "error at index: " << idx;
    }
}

/// Integration test to build a cylinder volume with contained surfaces
TEST(detector, detector_volume_construction) {

    using namespace detray;

    using detector_t =
        detray::detector<detector_registry::default_detector, covfie::field>;
    using transform3 = typename detector_t::transform3;
    using geo_obj_id = typename detector_t::geo_obj_ids;
    using mask_id = typename detector_t::masks::id;

    // portal factories
    using portal_cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_portal_cylinder2,
                        surface_id::e_portal>;
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

    // volume builder
    volume_builder<detector_t> vbuilder{};
    vbuilder.init_vol(d, volume_id::e_cylinder,
                      {0.f, 500.f, -1500.f, 1500.f, -M_PI, M_PI});
    const auto& vol = d.volumes().back();

    // initial checks
    EXPECT_EQ(d.volumes().size(), 2UL);
    EXPECT_EQ(d.surfaces().size(), 3UL);
    EXPECT_TRUE(vol.empty());
    EXPECT_EQ(vol.index(), 1UL);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);

    //
    // Fill the surface factories with data
    //

    // portal surfaces
    auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();
    typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
    // Creates two cylinders at radius 0mm and 10mm with an extent in z
    // of -5mm to 5mm and linking to volumes 0 and 2, respectively.
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), 0UL,
            std::vector<scalar>{10.f, -1500.f, 1500.f}));
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), 2UL,
            std::vector<scalar>{20.f, -1500.f, 1500.f}));
    pt_cyl_factory->add_components(std::move(cyl_sf_data));

    auto pt_disc_factory = std::make_shared<portal_disc_factory>();
    typename portal_disc_factory::sf_data_collection ring_sf_data;
    // Creates two discs with a radius of 10mm, linking to volumes 3 and 4
    ring_sf_data.push_back(std::make_unique<surface_data<detector_t, ring2D<>>>(
        transform3(point3{0.f, 0.f, -1500.f}), 3UL,
        std::vector<scalar>{0.f, 10.f}));
    ring_sf_data.push_back(std::make_unique<surface_data<detector_t, ring2D<>>>(
        transform3(point3{0.f, 0.f, 1500.f}), 4UL,
        std::vector<scalar>{0.f, 10.f}));
    pt_disc_factory->add_components(std::move(ring_sf_data));

    // sensitive surfaces
    auto ann_factory = std::make_shared<annulus_factory>();
    typename annulus_factory::sf_data_collection ann_sf_data;
    ann_sf_data.push_back(
        std::make_unique<surface_data<detector_t, annulus2D<>>>(
            transform3(point3{0.f, 0.f, 1000.f}), vol.index(),
            std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f}));
    ann_sf_data.push_back(
        std::make_unique<surface_data<detector_t, annulus2D<>>>(
            transform3(point3{0.f, 0.f, 1000.f}), vol.index(),
            std::vector<scalar>{350.f, 400.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f}));
    ann_factory->add_components(std::move(ann_sf_data));

    auto rect_factory = std::make_shared<rectangle_factory>();
    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, -10.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, -20.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, -30.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_factory->add_components(std::move(rect_sf_data));

    auto trpz_factory = std::make_shared<trapezoid_factory>();
    typename trapezoid_factory::sf_data_collection trpz_sf_data;
    trpz_sf_data.push_back(
        std::make_unique<surface_data<detector_t, trapezoid2D<>>>(
            transform3(point3{0.f, 0.f, 1000.f}), vol.index(),
            std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}));
    trpz_factory->add_components(std::move(trpz_sf_data));

    // passive surfaces
    auto cyl_factory = std::make_shared<cylinder_factory>();
    cyl_sf_data.clear();
    cyl_sf_data.push_back(
        std::make_unique<surface_data<detector_t, cylinder2D<>>>(
            transform3(point3{0.f, 0.f, 0.f}), vol.index(),
            std::vector<scalar>{5.f, -1300.f, 1300.f}));
    cyl_factory->add_components(std::move(cyl_sf_data));

    auto rng_factory = std::make_shared<ring_factory>();
    ring_sf_data.clear();
    ring_sf_data.push_back(std::make_unique<surface_data<detector_t, ring2D<>>>(
        transform3(point3{0.f, 0.f, -1300.f}), vol.index(),
        std::vector<scalar>{0.f, 5.f}));
    ring_sf_data.push_back(std::make_unique<surface_data<detector_t, ring2D<>>>(
        transform3(point3{0.f, 0.f, 1300.f}), vol.index(),
        std::vector<scalar>{0.f, 5.f}));
    rng_factory->add_components(std::move(ring_sf_data));

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
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, 10.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, 20.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_factory->add_components(std::move(rect_sf_data));
    rect_sf_data.clear();
    rect_sf_data.push_back(
        std::make_unique<surface_data<detector_t, rectangle2D<>>>(
            transform3(point3{0.f, 0.f, 30.f}), vol.index(),
            std::vector<scalar>{10.f, 8.f}));
    rect_factory->add_components(std::move(rect_sf_data));

    vbuilder.add_sensitives(rect_factory, geo_ctx);

    //
    // Adds all surfaces to the detector
    //
    vbuilder.build(d);

    //
    // check results
    //
    EXPECT_FALSE(vol.empty());
    // default detector makes no distinction between the surface types
    typename detector_registry::default_detector::object_link_type sf_range{};
    sf_range[0] = {3UL, 19UL};
    EXPECT_EQ(vol.full_range(), sf_range[geo_obj_id::e_sensitive]);
    EXPECT_EQ(vol.template obj_link<geo_obj_id::e_portal>(),
              sf_range[geo_obj_id::e_portal]);
    EXPECT_EQ(vol.template obj_link<geo_obj_id::e_sensitive>(),
              sf_range[geo_obj_id::e_sensitive]);
    EXPECT_EQ(vol.template obj_link<geo_obj_id::e_passive>(),
              sf_range[geo_obj_id::e_passive]);

    EXPECT_EQ(d.surfaces().size(), 19UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 3UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 4UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 3UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 3UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 7UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 4UL);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 2UL);

    // check surface type and volume link
    std::vector<surface_id> sf_ids{};
    sf_ids.reserve(d.surfaces().size());
    sf_ids.insert(sf_ids.end(), 3UL, surface_id::e_sensitive);
    sf_ids.insert(sf_ids.end(), 4UL, surface_id::e_portal);
    sf_ids.insert(sf_ids.end(), 6UL, surface_id::e_sensitive);
    sf_ids.insert(sf_ids.end(), 3UL, surface_id::e_passive);
    sf_ids.insert(sf_ids.end(), 3UL, surface_id::e_sensitive);
    std::vector<dindex> volume_links{};
    volume_links.insert(volume_links.end(), 3UL, 0UL);
    volume_links.insert(volume_links.end(), 16UL, 1UL);
    volume_links.reserve(d.surfaces().size());
    for (const auto [idx, sf_id] : detray::views::enumerate(sf_ids)) {
        const auto& sf = d.surface_by_index(idx);
        EXPECT_EQ(sf.id(), sf_id) << "error at index: " << idx;
        EXPECT_EQ(sf.volume(), volume_links.at(idx))
            << "error at index: " << idx;
    }

    // check that the transform indices are continuous
    for (std::size_t idx :
         detray::views::iota(d.transform_store().size() - 1)) {
        EXPECT_EQ(d.surface_by_index(idx).transform(), idx)
            << "error at index: " << idx;
    }

    // check surface mask links
    std::vector<typename detector_t::mask_link> mask_links{
        {mask_id::e_rectangle2, 0UL},
        {mask_id::e_annulus2, 0UL},
        {mask_id::e_trapezoid2, 0UL},
        {mask_id::e_portal_cylinder2, 0UL},
        {mask_id::e_portal_cylinder2, 1UL},
        {mask_id::e_portal_ring2, 0UL},
        {mask_id::e_portal_ring2, 1UL},
        {mask_id::e_annulus2, 1UL},
        {mask_id::e_annulus2, 2UL},
        {mask_id::e_rectangle2, 1UL},
        {mask_id::e_rectangle2, 2UL},
        {mask_id::e_rectangle2, 3UL},
        {mask_id::e_trapezoid2, 1UL},
        {mask_id::e_cylinder2, 2UL},
        {mask_id::e_ring2, 2UL},
        {mask_id::e_ring2, 3UL},
        {mask_id::e_rectangle2, 4UL},
        {mask_id::e_rectangle2, 5UL},
        {mask_id::e_rectangle2, 6UL}};
    for (const auto [idx, m_link] : detray::views::enumerate(mask_links)) {
        EXPECT_EQ(d.surface_by_index(idx).mask(), m_link)
            << "error at index: " << idx;
    }

    // check mask volume links
    volume_links.clear();
    volume_links = {0UL, 2UL, 1UL};
    check_mask<detector_t, mask_id::e_portal_cylinder2>(d, volume_links);
    check_mask<detector_t, mask_id::e_cylinder2>(d, volume_links);

    volume_links.clear();
    volume_links = {3UL, 4UL, 1UL, 1UL};
    check_mask<detector_t, mask_id::e_portal_ring2>(d, volume_links);
    check_mask<detector_t, mask_id::e_ring2>(d, volume_links);

    volume_links.clear();
    volume_links = {0UL, 1UL, 1UL};
    check_mask<detector_t, mask_id::e_annulus2>(d, volume_links);

    volume_links.clear();
    volume_links.reserve(7);
    volume_links.push_back(0UL);
    volume_links.insert(volume_links.end(), 6UL, 1UL);
    check_mask<detector_t, mask_id::e_rectangle2>(d, volume_links);

    volume_links.clear();
    volume_links = {0UL, 1UL};
    check_mask<detector_t, mask_id::e_trapezoid2>(d, volume_links);
}