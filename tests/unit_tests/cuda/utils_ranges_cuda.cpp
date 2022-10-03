/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "utils_ranges_cuda_kernel.hpp"
#include "vecmem/utils/cuda/copy.hpp"

using namespace detray;

// This tests the single value view
TEST(utils_ranges_cuda, single) {

    dindex value{251UL};
    dindex check{std::numeric_limits<dindex>::max()};

    // Run the test code
    single(value, check);

    // Check result value
    ASSERT_EQ(value, check);
}

// This tests the convenience range_enumeration function: range
/*TEST(utils_ranges_cuda, sequence_range) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input reference vector for test
    vecmem::vector<dindex> reference({2, 3, 4, 5, 6, 7});

    // Const range for enumeration test
    const darray<dindex, 2> range = {2, 7};

    // Output vector buffer for enumeration test
    vecmem::data::vector_buffer<dindex> check_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(
            range[1] - range[0] + 1),
        0, managed_resource);
    copy.setup(check_buffer);

    // Run test function
    sequence_range(range, check_buffer);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> check{&managed_resource};
    copy(check_buffer, check);

    // Check the result
    ASSERT_EQ(check, reference);
}

// This tests the convenience enumeration function
TEST(utils_ranges_cuda, enumerate_sequence) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<uint_holder> seq({{0}, {1}, {2}, {3}, {4}, {5}},
                                    &managed_resource);

    // Output vector buffer for enumeration test
    vecmem::data::vector_buffer<dindex> idx_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
        0, managed_resource);
    copy.setup(idx_buffer);

    vecmem::data::vector_buffer<unsigned int> uint_buffer(
        static_cast<vecmem::data::vector_buffer<unsigned int>::size_type>(
            seq.size()),
        0, managed_resource);
    copy.setup(uint_buffer);

    // Get vector_data object
    auto seq_data = vecmem::get_data(seq);

    // Run test function
    enumerate_sequence(idx_buffer, uint_buffer, seq_data);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> idx_vec{&managed_resource};
    copy(idx_buffer, idx_vec);

    vecmem::vector<unsigned int> uint_vec{&managed_resource};
    copy(uint_buffer, uint_vec);

    // Check the result
    for (unsigned int i = 0; i < idx_vec.size(); i++) {
        ASSERT_EQ(idx_vec[i], static_cast<dindex>(uint_vec[i]));
    }
}

// This tests the restricted iterator
TEST(utils_ranges_cuda, range) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<int> seq({0, 1, 2, 3, 4, 5}, &managed_resource);

    // Begin and end index for iteration
    const size_t begin = 1;
    const size_t end = 4;

    // Output vector buffer for iteration test
    vecmem::data::vector_buffer<int> check_buffer(
        static_cast<vecmem::data::vector_buffer<int>::size_type>(begin - end),
        0, managed_resource);
    copy.setup(check_buffer);

    // Get vector_data object
    auto seq_data = vecmem::get_data(seq);

    // Run test function
    iterate_range(check_buffer, seq_data, begin, end);

    // Copy vector buffer to output vector
    vecmem::vector<int> check{&managed_resource};
    copy(check_buffer, check);

    // Check the result
    ASSERT_EQ(check[0], 1);
    ASSERT_EQ(check[1], 2);
    ASSERT_EQ(check[2], 3);
}*/
