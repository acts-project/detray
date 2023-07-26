/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/tools/material_builder.hpp"
#include "detray/tools/material_factory.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>
#include <memory>
#include <vector>

using namespace detray;

namespace {

using point3 = __plugin::point3<scalar>;

using detector_t = detector<>;

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

}  // anonymous namespace

/// Unittest: Test the construction of a collection of materials
TEST(detray_tools, material_builder) {

    using transform3 = typename detector_t::transform3;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;

    // Build rectangle surfaces with material slabs
    using rectangle_factory =
        surface_factory<detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;
    auto mat_factory = std::make_unique<material_factory<detector_t>>(
        std::make_unique<rectangle_factory>());

    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    // auto geo_ctx = typename detector_t::geometry_context{};

    EXPECT_TRUE(d.material_store().template empty<material_id::e_slab>());
    EXPECT_TRUE(d.material_store().template empty<material_id::e_rod>());

    EXPECT_EQ(mat_factory->size(), 0u);
    EXPECT_TRUE(mat_factory->materials().empty());
    EXPECT_TRUE(mat_factory->thickness().empty());

    // Add material for a few rectangle surfaces
    mat_factory->push_back({transform3(point3{0.f, 0.f, -1.f}), 1u,
                            std::vector<scalar>{10.f, 8.f}});
    mat_factory->add_material(material_id::e_slab,
                              {1.f * unit<scalar>::mm, silicon<scalar>()});
    mat_factory->push_back({transform3(point3{0.f, 0.f, 1.f}), 1u,
                            std::vector<scalar>{20.f, 16.f}});
    mat_factory->add_material(material_id::e_slab,
                              {10.f * unit<scalar>::mm, tungsten<scalar>()});
    // Pass the parameters for 'gold'
    mat_factory->push_back({transform3(point3{0.f, 0.f, 1.f}), 1u,
                            std::vector<scalar>{20.f, 16.f}});
    mat_factory->add_material(
        material_id::e_rod,
        {0.1f * unit<scalar>::mm,
         std::vector<scalar>{
             3.344f * unit<scalar>::mm, 101.6f * unit<scalar>::mm, 196.97f, 79,
             19.32f * unit<scalar>::g / (1.f * unit<scalar>::cm3)},
         material_state::e_solid});

    EXPECT_EQ(mat_factory->size(), 3u);

    // Test the mayerial data
    EXPECT_NEAR(mat_factory->thickness()[0], 1.f * unit<scalar>::mm, tol);
    EXPECT_NEAR(mat_factory->thickness()[1], 10.f * unit<scalar>::mm, tol);
    EXPECT_NEAR(mat_factory->thickness()[2], 0.1f * unit<scalar>::mm, tol);
    EXPECT_EQ(mat_factory->materials()[0], silicon<scalar>());
    EXPECT_EQ(mat_factory->materials()[1], tungsten<scalar>());
    EXPECT_EQ(mat_factory->materials()[2], gold<scalar>());
}

/// Integration test: material builder as volume builder decorator
GTEST_TEST(detray_tools, decorator_material_builder) {

    using transform3 = typename detector_t::transform3;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;

    using portal_cylinder_factory_t =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_portal_cylinder2,
                        surface_id::e_portal>;
    using rectangle_factory =
        surface_factory<detector_t, rectangle2D<>, mask_id::e_rectangle2,
                        surface_id::e_sensitive>;
    using trapezoid_factory =
        surface_factory<detector_t, trapezoid2D<>, mask_id::e_trapezoid2,
                        surface_id::e_sensitive>;
    using cylinder_factory =
        surface_factory<detector_t, cylinder2D<>, mask_id::e_cylinder2,
                        surface_id::e_passive>;

    using mat_factory_t = material_factory<detector_t>;

    vecmem::host_memory_resource host_mr;
    detector_t d(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};

    auto vbuilder = std::make_unique<volume_builder<detector_t>>();
    auto mat_builder = material_builder<detector_t>{std::move(vbuilder)};

    EXPECT_TRUE(d.volumes().size() == 0);

    // Now init the volume
    mat_builder.init_vol(d, volume_id::e_cylinder);

    const auto &vol = d.volumes().back();
    EXPECT_TRUE(d.volumes().size() == 1u);
    EXPECT_EQ(vol.index(), 0u);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);

    // Add some portals first
    auto pt_cyl_factory = std::make_unique<portal_cylinder_factory_t>();

    pt_cyl_factory->push_back({transform3(point3{0.f, 0.f, 0.f}), 1u,
                               std::vector<scalar>{10.f, -1500.f, 1500.f}});
    pt_cyl_factory->push_back({transform3(point3{0.f, 0.f, 0.f}), 2u,
                               std::vector<scalar>{20.f, -1500.f, 1500.f}});

    // Then some passive and sensitive surfaces
    auto rect_factory = std::make_unique<rectangle_factory>();

    typename rectangle_factory::sf_data_collection rect_sf_data;
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -10.f}), 0u,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -20.f}), 0u,
                              std::vector<scalar>{10.f, 8.f});
    rect_sf_data.emplace_back(transform3(point3{0.f, 0.f, -30.f}), 0u,
                              std::vector<scalar>{10.f, 8.f});
    rect_factory->push_back(std::move(rect_sf_data));

    auto trpz_factory = std::make_unique<trapezoid_factory>();

    trpz_factory->push_back({transform3(point3{0.f, 0.f, 1000.f}), 0u,
                             std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}});

    auto cyl_factory = std::make_unique<cylinder_factory>();

    cyl_factory->push_back({transform3(point3{0.f, 0.f, 0.f}), 0u,
                            std::vector<scalar>{5.f, -1300.f, 1300.f}});

    // Now add the material for each surface
    auto mat_pt_cyl_factory =
        std::make_shared<mat_factory_t>(std::move(pt_cyl_factory));
    mat_pt_cyl_factory->add_material(
        material_id::e_slab, {1.f * unit<scalar>::mm, silicon<scalar>()});
    mat_pt_cyl_factory->add_material(
        material_id::e_slab, {1.5f * unit<scalar>::mm, silicon<scalar>()});

    auto mat_rect_factory =
        std::make_shared<mat_factory_t>(std::move(rect_factory));
    mat_rect_factory->add_material(material_id::e_slab,
                                   {1.f * unit<scalar>::mm, silicon<scalar>()});
    mat_rect_factory->add_material(material_id::e_slab,
                                   {2.f * unit<scalar>::mm, silicon<scalar>()});
    mat_rect_factory->add_material(material_id::e_slab,
                                   {3.f * unit<scalar>::mm, silicon<scalar>()});

    auto mat_trpz_factory =
        std::make_shared<mat_factory_t>(std::move(trpz_factory));
    mat_trpz_factory->add_material(material_id::e_slab,
                                   {1.f * unit<scalar>::mm, silicon<scalar>()});

    auto mat_cyl_factory =
        std::make_shared<mat_factory_t>(std::move(cyl_factory));
    mat_cyl_factory->add_material(
        material_id::e_slab, {1.5f * unit<scalar>::mm, tungsten<scalar>()});

    // Add surfaces and material to detector
    mat_builder.add_portals(mat_pt_cyl_factory, geo_ctx);
    mat_builder.add_sensitives(mat_rect_factory, geo_ctx);
    mat_builder.add_sensitives(mat_trpz_factory, geo_ctx);
    mat_builder.add_passives(mat_cyl_factory, geo_ctx);

    // Add the volume to the detector
    mat_builder.build(d);

    //
    // check results
    //
    EXPECT_EQ(d.surface_lookup().size(), 7u);
    EXPECT_EQ(d.transform_store().size(), 8u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(), 2u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 3u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1u);

    EXPECT_EQ(d.material_store().template size<material_id::e_slab>(), 7u);
    EXPECT_EQ(d.material_store().template size<material_id::e_rod>(), 0u);

    for (auto [idx, sf_desc] : detray::views::enumerate(d.surface_lookup())) {
        const auto &mat_link = sf_desc.material();
        EXPECT_EQ(mat_link.id(), material_id::e_slab);
        EXPECT_EQ(mat_link.index(), idx);
    }

    for (const auto &mat_slab :
         d.material_store().template get<material_id::e_slab>()) {
        EXPECT_TRUE(mat_slab.get_material() == silicon<scalar>() or
                    mat_slab.get_material() == tungsten<scalar>());
    }
}
