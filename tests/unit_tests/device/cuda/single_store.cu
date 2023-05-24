/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

#include "detray/core/detail/single_store.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using store_t = single_store<int, vecmem::vector>;
using store_dev_t = single_store<int, vecmem::device_vector>;

__global__ void basic_kernel(typename store_t::view_type view) {
    std::size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;

    store_dev_t store(view);

    if (globalIdx >= store.size()) {
        return;
    }

    printf("%d\n", store.at(globalIdx, empty_context{}));
}

GTEST_TEST(detray_detail, single_store) {

    // Run the kernel by copying the data to device

    vecmem::host_memory_resource host_mr;
    store_t store(host_mr);
    store.push_back(3);
    store.push_back(4);

    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy cpy;

    auto store_buffer = detray::get_buffer(store, dev_mr, cpy);

    basic_kernel<<<8, 1>>>(detray::get_data(store_buffer));

    // Run the same kernel with managed memory

    vecmem::cuda::managed_memory_resource mng_mr;
    store_t managed_store(mng_mr);
    managed_store.push_back(5);
    managed_store.push_back(6);

    basic_kernel<<<8, 1>>>(detray::get_data(managed_store));
}
