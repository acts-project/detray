/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Project include(s)
#include "ca_transform_store_cuda_kernel.hpp"
#include "detray/definitions/detail/algebra.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <climits>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

using namespace detray;

TEST(ca_transform_store_cuda, ca_transform_store) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    host_ca_transform_store_t static_ca_store(mng_mr);
    typename host_ca_transform_store_t::context_type ctx_0{0};
    typename host_ca_transform_store_t::context_type ctx_1{1};

    // Construct several transforms
    point3 t0{1.f, 2.f, 3.f};
    point3 z0{4.f, 2.f, 3.f};
    point3 x0{1.f, 0.f, 3.f};
    transform3 tf0{t0, z0, x0};

    point3 t1{1.f, 1.f, 2.f};
    point3 z1{1.f, 0.f, 0.f};
    point3 x1{0.f, 0.f, 2.f};
    transform3 tf1{t1, z1, x1};

    point3 t2{2.f, 2.f, 5.f};
    point3 z2{0.f, 0.f, 0.f};
    point3 x2{1.f, 2.f, 3.f};
    transform3 tf2{t2, z2, x2};

    point3 t3{2.f, 0.f, 5.f};
    point3 z3{1.f, 2.f, 3.f};
    point3 x3{0.f, 0.f, 0.f};
    transform3 tf3{t3, z3, x3};

    point3 t4{2.f, 0.f, 5.f};
    point3 z4{5.f, 3.f, 2.f};
    point3 x4{0.f, 1.f, 1.f};

    // Record the transforms into the store under two contexts
    // Testing push_back() and emplace_back()
    static_ca_store.reserve(3, ctx_0);
    static_ca_store.reserve(2, ctx_1);

    static_ca_store.push_back(tf0, ctx_0);
    static_ca_store.push_back(tf1, ctx_1);
    static_ca_store.push_back(tf2, ctx_0);
    static_ca_store.push_back(tf3, ctx_1);
    static_ca_store.emplace_back(ctx_0, t4, z4, x4);

    ASSERT_EQ(static_ca_store.size(ctx_0), 3u);
    ASSERT_EQ(static_ca_store.size(ctx_1), 2u);

    // Compare the transform operation results
    dindex n_transforms_ca0 = static_ca_store.size(ctx_0);
    dindex n_transforms_ca1 = static_ca_store.size(ctx_1);
    dindex n_transforms_total = n_transforms_ca0 + n_transforms_ca1;

    // Fill input vector
    vecmem::vector<point3> input(&mng_mr);
    for (unsigned int i = 0u; i < n_transforms_total; ++i) {
        input.push_back({static_cast<scalar>(i), static_cast<scalar>(i + 1u),
                         static_cast<scalar>(i + 2u)});
    }

    // Get transformed vectors from the host side
    std::vector<point3> output_host_ca0, output_host_ca1;
    auto range_ca0 = ranges::subrange(static_ca_store.get(ctx_0),
                                      dindex_range{0u, n_transforms_ca0});
    auto range_ca1 = ranges::subrange(static_ca_store.get(ctx_1),
                                      dindex_range{0u, n_transforms_ca1});

    std::size_t count{0u};
    for (const auto& tf : range_ca0) {
        auto output = tf.point_to_global(input[2 * count]);
        output_host_ca0.push_back(output);
        count++;
    }

    count = 0u;
    for (const auto& tf : range_ca1) {
        auto output = tf.point_to_global(input[1 + 2 * count]);
        output_host_ca1.push_back(output);
        count++;
    }

    // Get transformed vectors from the device side
    vecmem::vector<point3> output_device_ca0(n_transforms_ca0, &mng_mr);
    vecmem::vector<point3> output_device_ca1(n_transforms_ca1, &mng_mr);

    auto input_data = vecmem::get_data(input);
    auto static_ca_store_data = vecmem::get_data(*(static_ca_store.data()));

    // Context: 0
    auto output_data_ca0 = vecmem::get_data(output_device_ca0);
    ca_transform_test(input_data, static_ca_store_data, output_data_ca0,
                      static_ca_store.size(ctx_0), ctx_0);

    // Context: 1
    auto output_data_ca1 = vecmem::get_data(output_device_ca1);
    ca_transform_test(input_data, static_ca_store_data, output_data_ca1,
                      static_ca_store.size(ctx_1), ctx_1);

    // Compare the results of the transformations on host and device
    for (unsigned int i = 0u; i < n_transforms_ca0; ++i) {
        ASSERT_EQ(output_host_ca0[i], output_device_ca0[i]);
    }
    for (unsigned int i = 0u; i < n_transforms_ca1; ++i) {
        ASSERT_EQ(output_host_ca1[i], output_device_ca1[i]);
    }
}
