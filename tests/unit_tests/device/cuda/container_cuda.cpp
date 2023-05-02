/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "container_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <numeric>

using namespace detray;

TEST(container_cuda, multi_type_store) {

    // Vecmem memory resource
    vecmem::cuda::managed_memory_resource resource;

    // Create tuple vector store
    host_store_type store(resource);

    // Base store function check
    EXPECT_EQ(store.n_collections(), 3u);
    EXPECT_EQ(store.empty<0>(), true);
    EXPECT_EQ(store.empty<1>(), true);
    EXPECT_EQ(store.empty<2>(), true);

    // Add elements to the store
    store.emplace_back<0>(empty_context{}, 1UL);
    store.emplace_back<0>(empty_context{}, 2UL);
    store.emplace_back<1>(empty_context{}, 3.1f);
    store.emplace_back<1>(empty_context{}, 4.5f);
    store.emplace_back<2>(empty_context{}, 5.5);
    store.emplace_back<2>(empty_context{}, 6.0);

    vecmem::vector<std::size_t> int_vec{3UL, 4UL, 5UL};
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
    typename host_store_type::view_type store_data = get_data(store);

    vecmem::vector<double> cuda_sum(&resource);
    cuda_sum.push_back(0);
    auto sum_data = vecmem::get_data(cuda_sum);

    get_sum(store_data, sum_data);

    EXPECT_NEAR(cpu_sum, cuda_sum[0], 1e-6);
}