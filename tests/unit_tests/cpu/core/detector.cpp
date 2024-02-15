/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"

#include "detray/definitions/detail/indexing.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/test/utils/prefill_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

/// This tests the functionality of a detector as a data store manager
GTEST_TEST(detray_core, detector) {

    using namespace detray;

    using detector_t = detector<>;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using finder_id = typename detector_t::accel::id;

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
    EXPECT_TRUE(
        d.accelerator_store().template empty<finder_id::e_brute_force>());
    EXPECT_TRUE(d.accelerator_store().template empty<finder_id::e_disc_grid>());
    EXPECT_TRUE(
        d.accelerator_store().template empty<finder_id::e_cylinder2_grid>());
    EXPECT_TRUE(
        d.accelerator_store().template empty<finder_id::e_irr_disc_grid>());
    EXPECT_TRUE(d.accelerator_store()
                    .template empty<finder_id::e_irr_cylinder2_grid>());
    EXPECT_TRUE(d.accelerator_store().template empty<finder_id::e_default>());

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
    EXPECT_EQ(d.accelerator_store().template size<finder_id::e_brute_force>(),
              1u);
    EXPECT_EQ(d.accelerator_store().template size<finder_id::e_disc_grid>(),
              0u);
    EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_cylinder2_grid>(), 0u);
    EXPECT_EQ(d.accelerator_store().template size<finder_id::e_irr_disc_grid>(),
              0u);
    EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_irr_cylinder2_grid>(),
        0u);
    EXPECT_EQ(d.accelerator_store().template size<finder_id::e_default>(), 1u);
}
