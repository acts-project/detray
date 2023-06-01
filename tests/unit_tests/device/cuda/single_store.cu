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
#include "detray/core/detail/tuple_container copy.hpp"
#include "detray/definitions/cuda_definitions.hpp"

// GTest include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

constexpr float tol{1e-6f};

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


using t0 = vecmem_types<unsigned, true>; // vector of unsigned
using t1 = vecmem_types<int, false>; // 1 int
using t2 = vecmem_types<float, true>; // vector of float
using t3 = vecmem_types<double, false>; // 1 double

using test_tuple_ts = detray::detail::tuple_types<t0,t1,t2,t3>;


__global__ void basic_kernel(test_tuple_ts::view tuple_view) {
    std::size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;

    if(globalIdx != 0) {
        return;
    }
    
    test_tuple_ts::device tuple_device(tuple_view);

    auto uns_vec =   get<0>(tuple_device);
    auto& _int =      *get<1>(tuple_device);
    auto float_vec = get<2>(tuple_device);
    auto& _double =   *get<3>(tuple_device);

    uns_vec.at(0) += 1;
    uns_vec.at(1) += 3;
    _int += 1;
    float_vec.at(0) += 1.1;
    float_vec.at(1) += 1.1;
    _double += 2.2;
}

GTEST_TEST(detray_detail, tuple_types) {
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy copy;

    test_tuple_ts::host host(vecmem::vector<unsigned>(1, 0), 1, vecmem::vector<float>(2, 2.1), 3.3);

    EXPECT_EQ(detail::get<0>(host).size(), 1);
    EXPECT_EQ(detail::get<0>(host)[0], 0);
    EXPECT_EQ(detail::get<1>(host), 1);
    EXPECT_EQ(detail::get<2>(host).size(), 2);
    EXPECT_NEAR(detail::get<2>(host)[0], 2.1, tol);
    EXPECT_NEAR(detail::get<2>(host)[1], 2.1, tol);
    EXPECT_EQ(detail::get<3>(host), 3.3);

    // Edit data on host type
    detail::get<0>(host).push_back(1);
    EXPECT_EQ(detail::get<0>(host).size(), 2);
    EXPECT_EQ(detail::get<0>(host)[1], 1);

    detail::get<1>(host) -= 2;
    EXPECT_EQ(detail::get<1>(host), -1);
    
    // Can also create view and device types on top of a host type
    test_tuple_ts::view host_view(host);
    test_tuple_ts::device host_device(host_view);
    
    // And edit using the device type
    detail::get<2>(host_device)[1] += 0.1;
    EXPECT_NEAR(detail::get<2>(host)[1], 2.2, tol);

    // Create buffer on device
    test_tuple_ts::buffer dev_buff(dev_mr, detail::get<0>(host).size(), detail::get<2>(host).size());
    test_tuple_ts::view dev_view(dev_buff);

    // Copy the vecmem vectors
    copy(detail::get<0>(host_view), detail::get<0>(dev_view));
    copy(detail::get<2>(host_view), detail::get<2>(dev_view));
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // Copy PODs (note that in cudaMemcpy api src and target are inverse to vecmem api :/)
    cudaMemcpy(detail::get<1>(dev_view), detail::get<1>(host_view), sizeof(t1::value_type), cudaMemcpyHostToDevice);
    cudaMemcpy(detail::get<3>(dev_view), detail::get<3>(host_view), sizeof(t3::value_type), cudaMemcpyHostToDevice);
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // Edit data on device
    basic_kernel<<<1,1>>>(dev_view);
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // Copy results back to host
    test_tuple_ts::host endCopy;
    copy(detail::get<0>(dev_view), detail::get<0>(host_view));
    copy(detail::get<2>(dev_view), detail::get<2>(host_view));
    cudaMemcpy(detail::get<1>(host_view), detail::get<1>(dev_view), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(detail::get<3>(host_view), detail::get<3>(dev_view), sizeof(double), cudaMemcpyDeviceToHost);
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());


    EXPECT_EQ(detail::get<0>(host).size(), 2);
    EXPECT_EQ(detail::get<0>(host)[0], 1);
    EXPECT_EQ(detail::get<0>(host)[1], 4);
    EXPECT_EQ(detail::get<1>(host), 0);
    EXPECT_EQ(detail::get<2>(host).size(), 2);
    EXPECT_NEAR(detail::get<2>(host)[0], 3.2, tol);
    EXPECT_NEAR(detail::get<2>(host)[1], 3.3, tol);
    EXPECT_EQ(detail::get<3>(host), 5.5);
}