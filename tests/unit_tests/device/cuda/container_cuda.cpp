/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "container_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <iostream>
#include <numeric>

using namespace detray;

constexpr double tol{1e-6};

/// Test the access to the single store by summing the contained values
TEST(container_cuda, single_store) {

    // Vecmem memory resources
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy cpy;

    // Create single store(s)
    single_store_t store(host_mr);
    single_store_t mng_store(mng_mr);

    // Base store function check
    EXPECT_TRUE(store.empty());
    EXPECT_EQ(store.size(), 0u);
    EXPECT_TRUE(mng_store.empty());
    EXPECT_EQ(mng_store.size(), 0u);

    // Test the managed memory allocation
    empty_context ctx{};
    mng_store.reserve(4, ctx);
    mng_store.emplace_back(ctx, 1.);
    mng_store.push_back(2., ctx);
    mng_store.insert(vecmem::vector<double>{10.5, 7.6}, ctx);
    store.append(mng_store, ctx);

    EXPECT_FALSE(store.empty());
    EXPECT_EQ(store.size(), 4u);
    EXPECT_FALSE(mng_store.empty());
    EXPECT_EQ(mng_store.size(), 4u);

    // Check the host-side access to the data
    EXPECT_NEAR(mng_store[0], 1., tol);
    EXPECT_NEAR(mng_store[2], 10.5, tol);
    EXPECT_NEAR(mng_store.at(1, ctx), 2., tol);
    EXPECT_NEAR(mng_store.at(3, ctx), 7.6, tol);

    // CPU sum check
    double cpu_sum{std::accumulate(mng_store.begin(), mng_store.end(), 0.)};
    EXPECT_NEAR(cpu_sum, 21.1, tol);
    EXPECT_NEAR(cpu_sum, std::accumulate(store.begin(), store.end(), 0.), tol);

    // Get the view to the managed memory and run the test
    single_store_t::view_type mng_store_view = detray::get_data(mng_store);

    vecmem::vector<double> cuda_sum(&mng_mr);
    cuda_sum.push_back(0.);
    vecmem::data::vector_view<double> sum_data = vecmem::get_data(cuda_sum);

    test_single_store(mng_store_view, sum_data);

    EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);

    // Check that the store can be cleared
    mng_store.clear(ctx);

    EXPECT_TRUE(mng_store.empty());
    EXPECT_EQ(mng_store.size(), 0u);

    // Reset for next test
    cuda_sum[0] = 0.;

    // Copy the host store to device, get the view to it and run the test again
    single_store_t::buffer_type store_buffer =
        detray::get_buffer(store, dev_mr, cpy);
    single_store_t::view_type buffer_view = detray::get_data(store_buffer);

    test_single_store(buffer_view, sum_data);

    EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);
}

/*TEST(container_cuda, regular_multi_store) {

    // Vecmem memory resource
    vecmem::cuda::managed_memory_resource resource;

    // Create tuple vector store
    reg_multi_store_t store(resource);

    // Base store function check
    EXPECT_EQ(store.n_collections(), 3u);
    EXPECT_EQ(store.empty<0>(), true);
    EXPECT_EQ(store.empty<1>(), true);
    EXPECT_EQ(store.empty<2>(), true);

    // Add elements to the store
    store.emplace_back<0>(empty_context{}, 1u);
    store.emplace_back<0>(empty_context{}, 2u);
    store.emplace_back<1>(empty_context{}, 3.1f);
    store.emplace_back<1>(empty_context{}, 4.5f);
    store.emplace_back<2>(empty_context{}, 5.5);
    store.emplace_back<2>(empty_context{}, 6.0);

    vecmem::vector<std::size_t> int_vec{3u, 4u, 5u};
    store.insert(int_vec);

    vecmem::vector<float> float_vec{12.1f, 5.6f};
    store.insert(float_vec);

    store.insert(vecmem::vector<double>{10.5, 7.6});

    // CPU sum check
    double cpu_sum = 0.;
    cpu_sum =
        std::accumulate(store.get<0>().begin(), store.get<0>().end(), cpu_sum);
    cpu_sum =
        std::accumulate(store.get<1>().begin(), store.get<1>().end(), cpu_sum);
    cpu_sum =
        std::accumulate(store.get<2>().begin(), store.get<2>().end(), cpu_sum);
    EXPECT_NEAR(cpu_sum, 69.9, 1e-6);

    // CUDA sum check
    typename reg_multi_store_t::view_type store_data = get_data(store);

    vecmem::vector<double> cuda_sum(&resource);
    cuda_sum.push_back(0);
    auto sum_data = vecmem::get_data(cuda_sum);

    test_reg_multi_store(store_data, sum_data);

    EXPECT_NEAR(cpu_sum, cuda_sum[0], 1e-6);
}*/
