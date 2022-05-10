/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/mask_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/masks/masks.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a static transform store
TEST(ALGEBRA_PLUGIN, static_mask_store) {

    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    enum mask_ids : unsigned int {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_ring2 = 4,
        e_single3 = 5,
    };

    using rectangle = rectangle2<>;
    using trapezoid = trapezoid2<>;
    using annulus = annulus2<>;
    using cylinder = cylinder3<>;
    using ring = ring2<>;
    using single = single3<0>;

    // Types must be sorted according to their id (here: masks/mask_identifier)
    using mask_defs = mask_registry<mask_ids, rectangle, trapezoid, annulus,
                                    cylinder, ring, single>;
    using id = typename mask_defs::id;
    using mask_container_t = typename mask_defs::mask_store_type<>;
    mask_container_t store(host_mr);

    ASSERT_TRUE(store.empty<id::e_annulus2>());
    ASSERT_TRUE(store.empty<id::e_cylinder3>());
    ASSERT_TRUE(store.empty<id::e_rectangle2>());
    ASSERT_TRUE(store.empty<id::e_ring2>());
    ASSERT_TRUE(store.empty<id::e_single3>());
    ASSERT_TRUE(store.empty<id::e_trapezoid2>());

    store.add_value<id::e_cylinder3>(1., 0.5, 2.0, 0u);

    ASSERT_TRUE(store.empty<id::e_annulus2>());
    ASSERT_EQ(store.size<id::e_cylinder3>(), 1);
    ASSERT_TRUE(store.empty<id::e_rectangle2>());
    ASSERT_TRUE(store.empty<id::e_ring2>());
    ASSERT_TRUE(store.empty<id::e_single3>());
    ASSERT_TRUE(store.empty<id::e_trapezoid2>());

    store.add_value<id::e_cylinder3>(1., 1.5, 2.0, 0u);
    store.add_value<id::e_trapezoid2>(0.5, 1.5, 4.0, 0u);
    store.add_value<id::e_rectangle2>(1.0, 2.0, 0u);
    store.add_value<id::e_rectangle2>(2.0, 1.0, 0u);
    store.add_value<id::e_rectangle2>(10.0, 100.0, 0u);

    ASSERT_TRUE(store.empty<id::e_annulus2>());
    ASSERT_EQ(store.size<id::e_cylinder3>(), 2);
    ASSERT_EQ(store.size<id::e_rectangle2>(), 3);
    ASSERT_TRUE(store.empty<id::e_ring2>());
    ASSERT_TRUE(store.empty<id::e_single3>());
    ASSERT_EQ(store.size<id::e_trapezoid2>(), 1);

    store.add_value<id::e_annulus2>(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0u);
    store.add_value<id::e_ring2>(10.0, 100.0, 0u);
    store.add_value<id::e_ring2>(10.0, 100.0, 0u);
    store.add_value<id::e_ring2>(10.0, 100.0, 0u);
    store.add_value<id::e_ring2>(10.0, 100.0, 0u);

    const auto &annulus_masks = store.group<id::e_annulus2>();
    const auto &cylinder_masks = store.group<id::e_cylinder3>();
    const auto &rectangle_masks = store.group<id::e_rectangle2>();
    const auto &ring_masks = store.group<id::e_ring2>();
    const auto &single_masks = store.group<id::e_single3>();
    const auto &trapezoid_masks = store.group<id::e_trapezoid2>();

    ASSERT_TRUE(annulus_masks.size() == 1);
    ASSERT_TRUE(cylinder_masks.size() == 2);
    ASSERT_TRUE(rectangle_masks.size() == 3);
    ASSERT_TRUE(ring_masks.size() == 4);
    ASSERT_TRUE(single_masks.size() == 0);
    ASSERT_TRUE(trapezoid_masks.size() == 1);
}
