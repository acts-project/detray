/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "mask_store_cuda_kernel.hpp"

TEST(mask_store_cuda, mask_store) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Types must be sorted according to their id (here: masks/mask_identifier)
    mask_store<dtuple, vecmem::vector, rectangle, trapezoid, ring, cylinder,
               single, annulus>
        store(mng_mr);

    ASSERT_TRUE(store.template empty<e_annulus2>());
    ASSERT_TRUE(store.template empty<e_cylinder3>());
    ASSERT_TRUE(store.template empty<e_rectangle2>());
    ASSERT_TRUE(store.template empty<e_ring2>());
    ASSERT_TRUE(store.template empty<e_single3>());
    ASSERT_TRUE(store.template empty<e_trapezoid2>());

    store.template add_mask<e_cylinder3>(1., 0.5, 2.0);

    ASSERT_TRUE(store.template empty<e_annulus2>());
    ASSERT_EQ(store.template size<e_cylinder3>(), 1);
    ASSERT_TRUE(store.template empty<e_rectangle2>());
    ASSERT_TRUE(store.template empty<e_ring2>());
    ASSERT_TRUE(store.template empty<e_single3>());
    ASSERT_TRUE(store.template empty<e_trapezoid2>());

    store.template add_mask<e_cylinder3>(1., 1.5, 2.0);
    store.template add_mask<e_trapezoid2>(0.5, 1.5, 4.0);
    store.template add_mask<e_rectangle2>(1.0, 2.0);
    store.template add_mask<e_rectangle2>(2.0, 1.0);
    store.template add_mask<e_rectangle2>(10.0, 100.0);

    ASSERT_TRUE(store.template empty<e_annulus2>());
    ASSERT_EQ(store.template size<e_cylinder3>(), 2);
    ASSERT_EQ(store.template size<e_rectangle2>(), 3);
    ASSERT_TRUE(store.template empty<e_ring2>());
    ASSERT_TRUE(store.template empty<e_single3>());
    ASSERT_EQ(store.template size<e_trapezoid2>(), 1);

    store.template add_mask<e_annulus2>(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    store.template add_mask<e_ring2>(10.0, 100.0);
    store.template add_mask<e_ring2>(10.0, 100.0);
    store.template add_mask<e_ring2>(10.0, 100.0);
    store.template add_mask<e_ring2>(10.0, 100.0);

    const auto &annulus_masks = store.template group<e_annulus2>();
    const auto &cylinder_masks = store.template group<e_cylinder3>();
    const auto &rectangle_masks = store.template group<e_rectangle2>();
    const auto &ring_masks = store.template group<e_ring2>();
    const auto &single_masks = store.template group<e_single3>();
    const auto &trapezoid_masks = store.template group<e_trapezoid2>();

    ASSERT_TRUE(annulus_masks.size() == 1);
    ASSERT_TRUE(cylinder_masks.size() == 2);
    ASSERT_TRUE(rectangle_masks.size() == 3);
    ASSERT_TRUE(ring_masks.size() == 4);
    ASSERT_TRUE(single_masks.size() == 0);
    ASSERT_TRUE(trapezoid_masks.size() == 1);

    auto store_data = get_data(store);
    mask_test(store_data);
}
