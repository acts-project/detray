/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/masks/rectangle2.hpp"
#include "detray/materials/material_slab.hpp"
#include "tests/common/tools/create_single_plane_detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using vector3 = __plugin::point3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;
using registry_type = detector_registry::default_detector;

TEST(single_plane_detector, rectangle) {
    vecmem::host_memory_resource host_resource;

    constexpr registry_type::mask_ids mask_id =
        registry_type::mask_ids::e_rectangle2;
    constexpr registry_type::material_ids material_id =
        registry_type::material_ids::e_slab;

    const vector3 x{1, 0, 0};
    const vector3 z{0, 0, 1};
    const vector3 t{0, 0, 0};

    const scalar hx = 100. * unit_constants::mm;
    const scalar hy = 100. * unit_constants::mm;

    const material<scalar> mat = silicon<scalar>();
    const scalar thickness = 2 * unit_constants::mm;

    single_plane_detector_creator<mask_id, material_id> det_creator(
        host_resource);

    det_creator.set_transform(t, z, x);
    det_creator.set_mask(hx, hy);
    det_creator.set_material(mat, thickness);
    auto det = det_creator.build();

    auto masks = det.mask_store().group<mask_id>();
    auto materials = det.material_store().group<material_id>();

    ASSERT_EQ(masks.size(), 1);
    ASSERT_EQ(materials.size(), 1);
    EXPECT_FLOAT_EQ(masks[0][0], hx);
    EXPECT_FLOAT_EQ(masks[0][1], hy);

    EXPECT_EQ(materials[0].get_material(), mat);
    EXPECT_FLOAT_EQ(materials[0].thickness(), thickness);
}