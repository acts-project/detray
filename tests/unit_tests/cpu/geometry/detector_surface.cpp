/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/test/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

namespace {

/// Define mask types
enum mask_ids : unsigned int {
    e_unmasked = 0u,
};

/// Define material types
enum material_ids : unsigned int {
    e_slab = 0u,
};

constexpr detray::scalar tol{1e-6f};

}  // anonymous namespace

// This tests the construction of a surface descriptor object
GTEST_TEST(detray_geometry, surface_descriptor) {

    using namespace detray;

    using mask_link_t = dtyped_index<mask_ids, dindex>;
    using material_link_t = dtyped_index<material_ids, dindex>;

    mask_link_t mask_id{mask_ids::e_unmasked, 0u};
    material_link_t material_id{material_ids::e_slab, 0u};

    surface_descriptor<mask_link_t, material_link_t> desc(
        1u, mask_id, material_id, 2u, surface_id::e_sensitive);

    // Test access
    ASSERT_EQ(desc.transform(), 1u);
    ASSERT_EQ(desc.volume(), 2u);
    ASSERT_EQ(desc.id(), surface_id::e_sensitive);
    ASSERT_FALSE(desc.is_portal());
    ASSERT_FALSE(desc.is_passive());
    ASSERT_TRUE(desc.is_sensitive());
    ASSERT_EQ(desc.mask(), mask_id);
    ASSERT_EQ(desc.material(), material_id);

    // Test setters
    desc.set_volume(5u);
    desc.set_id(surface_id::e_portal);
    desc.set_index(6u);
    desc.update_transform(7u);
    desc.update_mask(7u);
    desc.update_material(7u);

    ASSERT_EQ(desc.transform(), 8u);
    ASSERT_EQ(desc.volume(), 5u);
    ASSERT_EQ(desc.id(), surface_id::e_portal);
    ASSERT_EQ(desc.index(), 6u);
    ASSERT_TRUE(desc.is_portal());
    ASSERT_FALSE(desc.is_passive());
    ASSERT_FALSE(desc.is_sensitive());
    ASSERT_EQ(desc.mask().id(), mask_ids::e_unmasked);
    ASSERT_EQ(desc.mask().index(), 7u);
    ASSERT_EQ(desc.material().id(), material_ids::e_slab);
    ASSERT_EQ(desc.material().index(), 7u);
}

/// This tests the functionality of a surface interface
GTEST_TEST(detray_geometry, surface) {

    using namespace detray;

    using detector_t = detector<toy_metadata>;

    using point2_t = surface<detector_t>::point2;
    using point3_t = surface<detector_t>::point3;
    using vector3_t = surface<detector_t>::vector3;

    vecmem::host_memory_resource host_mr;
    const auto [toy_det, names] = build_toy_detector(host_mr);

    auto ctx = typename detector_t::geometry_context{};

    const auto disc_descr = toy_det.surfaces()[1u];
    const auto disc = surface{toy_det, disc_descr};

    // IDs
    ASSERT_EQ(disc.barcode(), disc_descr.barcode());
    ASSERT_EQ(disc.volume(), 0u);
    ASSERT_EQ(disc.index(), 1u);
    ASSERT_EQ(disc.id(), surface_id::e_portal);
    ASSERT_EQ(disc.shape_id(), detector_t::masks::id::e_portal_ring2);
    ASSERT_FALSE(disc.is_sensitive());
    ASSERT_FALSE(disc.is_passive());
    ASSERT_TRUE(disc.is_portal());

    // Transformation matrix
    const auto disc_translation =
        vector3_t{0.f, 0.f, -824.5f};  // beampipe portal
    ASSERT_EQ(disc.transform(ctx).translation(), disc_translation);
    ASSERT_EQ(disc.center(ctx), disc_translation);

    // Surface normal
    const auto z_axis = vector3_t{0.f, 0.f, 1.f};
    // trigger all code paths
    ASSERT_EQ(disc.normal(ctx, point3_t{0.f, 0.f, 0.f}), z_axis);
    ASSERT_EQ(disc.normal(ctx, point2_t{0.f, 0.f}), z_axis);

    // Incidence angle
    const auto dir = vector3_t{1.f, 1.f, 1.f};
    ASSERT_NEAR(disc.cos_angle(ctx, dir, point3_t{0.f, 0.f, 0.f}), 1.f, tol);
    ASSERT_NEAR(disc.cos_angle(ctx, dir, point2_t{0.f, 0.f}), 1.f, tol);

    // Coordinate transformations
    const point3_t glob_pos = {4.f, 7.f, 4.f};
    const point3_t local = disc.global_to_local(ctx, glob_pos, {});
    const point2_t bound = disc.global_to_bound(ctx, glob_pos, {});

    ASSERT_NEAR(local[0], std::sqrt(65.f), tol);
    ASSERT_NEAR(local[1], std::atan2(7.f, 4.f), tol);
    ASSERT_NEAR(bound[0], local[0], tol);
    ASSERT_NEAR(bound[1], local[1], tol);

    // Roundtrip
    const point3_t global = disc.local_to_global(ctx, local, {});
    const point3_t global2 = disc.bound_to_global(ctx, bound, {});

    ASSERT_NEAR(glob_pos[0], global[0], tol);
    ASSERT_NEAR(glob_pos[1], global[1], tol);
    ASSERT_NEAR(glob_pos[2], global[2], tol);

    ASSERT_NEAR(global2[0], glob_pos[0], tol);
    ASSERT_NEAR(global2[1], glob_pos[1], tol);
    // The bound transform assumes the point is on surface
    ASSERT_NEAR(global2[2], disc_translation[2], tol);
}
