/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "tuple_test_cuda_kernel.hpp"

TEST(tuple_test_cuda, tuple_test) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // mask_store-like tuple container
    vec_tuple<thrust::tuple, dvector, int, float, double> input_host(mng_mr);
    detail::get<0>(input_host._tuple) = vecmem::vector<int>({1, 2, 3});
    detail::get<1>(input_host._tuple) = vecmem::vector<float>({1.1, 5, 6});
    detail::get<2>(input_host._tuple) = vecmem::vector<double>({2.1, 2.2, 0.});

    // Get vec_tuple_data with vecmem::data::vector_view
    auto input_data = get_data(input_host);

    // output vector which will copy the contents of input_host
    vecmem::vector<int> output1({0, 0, 0}, &mng_mr);
    vecmem::vector<float> output2({0, 0, 0}, &mng_mr);
    vecmem::vector<double> output3({0, 0, 0}, &mng_mr);
    auto output1_data = vecmem::get_data(output1);
    auto output2_data = vecmem::get_data(output2);
    auto output3_data = vecmem::get_data(output3);

    // run the cuda kernel for copy the tuple data into output vectors
    tuple_copy(input_data, output1_data, output2_data, output3_data);

    // Check if the copied vector is same with the input data
    ASSERT_EQ(detail::get<0>(input_host._tuple), output1);
    ASSERT_EQ(detail::get<1>(input_host._tuple), output2);
    ASSERT_EQ(detail::get<2>(input_host._tuple), output3);
}