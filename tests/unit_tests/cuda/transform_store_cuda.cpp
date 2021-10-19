/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "transform_store_cuda_kernel.hpp"

TEST(transform_store_cuda, transform_store) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    static_transform_store<> static_store;
    static_transform_store<>::context ctx0;
    static_transform_store<>::context ctx1;

    ASSERT_TRUE(static_store.empty(ctx0));

    ASSERT_EQ(static_store.size(ctx0), 0u);

    point3 t0{0., 0., 0.};
    transform3 tf0{t0};
    static_store.push_back(ctx0, tf0);
    ASSERT_EQ(static_store.size(ctx0), 1u);

    point3 t1{1., 0., 0.};
    transform3 tf1{t1};
    static_store.push_back(ctx1, tf1);
    ASSERT_EQ(static_store.size(ctx1), 2u);

    point3 t2{2., 0., 0.};
    transform3 tf2{t2};
    static_store.push_back(ctx0, std::move(tf2));
    ASSERT_EQ(static_store.size(ctx0), 3u);

    point3 t3{2., 0., 0.};
    static_store.emplace_back(ctx0, std::move(t3));
    ASSERT_EQ(static_store.size(ctx0), 4u);

    static_store.emplace_back(ctx0);
    ASSERT_EQ(static_store.size(ctx0), 5u);

    auto static_store_data = get_data(static_store);
    replace_test(static_store_data);
}
