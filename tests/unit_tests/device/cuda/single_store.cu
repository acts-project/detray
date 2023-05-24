/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

#include "detray/core/detail/single_store copy.hpp"

// GTest include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

using types = single_store_types<int>;

__global__ void basic_kernel(types::view view) {
    std::size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;

    if (globalIdx >= view.size()) {
        return;
    }

    view.at(globalIdx) += globalIdx;
}

GTEST_TEST(detray_detail, single_store) {

    types::host host(3, 3);
    host.at(1) = 100;
    host.push_back(-1);

    types::view host_view(host.get_data());
    host_view.at(2) -= 1;

    EXPECT_EQ(host.size(), 4);
    EXPECT_EQ(host.at(0), 3);
    EXPECT_EQ(host.at(1), 100);
    EXPECT_EQ(host.at(2), 2);
    EXPECT_EQ(host.at(3), -1);

    vecmem::host_memory_resource host_mr;
    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy cpy;

    types::buffer dev_buffer(host.get_data(), dev_mr, cpy);

    types::view dev_view(dev_buffer.get_data());

    basic_kernel<<<8, 1>>>(dev_view);

    types::buffer host_buffer(dev_buffer.get_data(), host_mr, cpy);
    types::view host_view_2(host_buffer.get_data());

    EXPECT_EQ(host_view_2.size(), 4);
    EXPECT_EQ(host_view_2.at(0), 3);
    EXPECT_EQ(host_view_2.at(1), 101);
    EXPECT_EQ(host_view_2.at(2), 4);
    EXPECT_EQ(host_view_2.at(3), 2);
}
