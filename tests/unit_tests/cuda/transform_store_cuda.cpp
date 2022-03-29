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

#include "transform_store_cuda_kernel.hpp"

TEST(transform_store_cuda, transform_store) {

    using point3 = __plugin::point3<detray::scalar>;
    using transform3 = __plugin::transform3<detray::scalar>;

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    static_transform_store<> static_store(mng_mr);
    static_transform_store<>::context ctx0{};
    static_transform_store<>::context ctx1{};

    ASSERT_TRUE(static_store.empty(ctx0));

    ASSERT_EQ(static_store.size(ctx0), 0u);

    point3 t0{1., 2., 3.};
    point3 z0{4., 2., 3.};
    point3 x0{1., 0, 3.};
    transform3 tf0{t0, z0, x0};
    static_store.push_back(ctx0, tf0);
    ASSERT_EQ(static_store.size(ctx0), 1u);

    point3 t1{1., 1., 2.};
    point3 z1{1., 0, 0};
    point3 x1{0, 0, 2.};
    transform3 tf1{t1, z1, x1};
    static_store.push_back(ctx1, tf1);
    ASSERT_EQ(static_store.size(ctx1), 2u);

    point3 t2{2., 2., 5.};
    point3 z2{0., 0., 0.};
    point3 x2{1., 2., 3.};
    transform3 tf2{t2, z2, x2};
    static_store.push_back(ctx0, std::move(tf2));
    ASSERT_EQ(static_store.size(ctx0), 3u);

    point3 t3{2., 0., 5.};
    point3 z3{1., 2., 3.};
    point3 x3{0., 0., 0.};
    transform3 tf3{t3, z3, x3};
    static_store.emplace_back(ctx0, std::move(t3));
    ASSERT_EQ(static_store.size(ctx0), 4u);

    static_store.emplace_back(ctx0);
    ASSERT_EQ(static_store.size(ctx0), 5u);

    // Compare the transform operation results
    unsigned int n_transforms = static_store.size(ctx0);

    // Fill input vector
    vecmem::vector<point3> input(&mng_mr);
    for (unsigned int i = 0; i < n_transforms; i++) {
        input.push_back({static_cast<scalar>(i), static_cast<scalar>(i + 1),
                         static_cast<scalar>(i + 2)});
    }

    // Get transformed vector from host side
    std::vector<point3> output_host;

    auto range = static_store.range(0, n_transforms, ctx0);
    int count = 0;
    for (auto tf : range) {
        auto output = tf.point_to_global(input[count]);
        output_host.push_back(output);
        count++;
    }

    // Get transformed vector from device side
    vecmem::vector<point3> output_device(n_transforms, &mng_mr);

    auto input_data = vecmem::get_data(input);
    auto output_data = vecmem::get_data(output_device);
    auto static_store_data = get_data(static_store);
    transform_test(input_data, static_store_data, output_data);

    // Compare the values
    for (unsigned int i = 0; i < n_transforms; i++) {
        ASSERT_EQ(output_host[i], output_device[i]);
    }
}
