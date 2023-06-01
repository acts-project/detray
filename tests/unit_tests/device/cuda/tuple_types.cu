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

#include "detray/core/detail/tuple_container copy.hpp"
#include "detray/definitions/cuda_definitions.hpp"

// GTest include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

constexpr float tol{1e-6f};

using t0 = vecmem_types<unsigned, true>; // vector of unsigned
using t1 = vecmem_types<int, false>; // 1 int
using t2 = vecmem_types<float, true>; // vector of float
using t3 = vecmem_types<double, false>; // 1 double

using test_tuple_ts = detray::detail::tuple_types<t0,t1,t2,t3>;


template<typename vec1_t, typename vec2_t>
__host__ __device__ float dosomethingcomplicated(vec1_t vec1, vec2_t vec2) {
    float result = 0.;
    for(std::size_t i = 0; i < vec1.size(); ++i) {
        result += static_cast<float>(vec1[i]);
    }
    for(std::size_t j = 0; j < vec2.size(); ++j) {
        result += static_cast<float>(vec2[j]);
    }
    return result;
}

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
    _double += dosomethingcomplicated(uns_vec, float_vec);
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
    EXPECT_NEAR(detail::get<3>(host), 3.3, tol);

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
    EXPECT_NEAR(detail::get<3>(host), 17., tol);

    detail::get<3>(host) -= dosomethingcomplicated(detail::get<0>(host), detail::get<2>(host));
    EXPECT_NEAR(detail::get<3>(host), 5.5, tol);

}