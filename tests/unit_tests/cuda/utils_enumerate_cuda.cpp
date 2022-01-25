/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>
//#include <vecmem/memory/cuda/device_memory_resource.hpp>
//#include <vecmem/memory/cuda/host_memory_resource.hpp>

#include "utils_enumerate_cuda_kernel.hpp"
#include "vecmem/utils/cuda/copy.hpp"

using namespace detray;

// This tests the convenience range_enumeration function: single
TEST(utils_enumerate_cuda, sequence_single) {

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    vecmem::vector<dindex> check(&managed_resource);
    check.push_back(0);

    vecmem::vector<dindex> single(&managed_resource);
    single.push_back(7);

    auto check_data = vecmem::get_data(check);
    auto single_data = vecmem::get_data(single);

    sequence_single(check_data, single_data);

    ASSERT_EQ(check[0], single[0]);
}

// This tests the convenience range_enumeration function: range
TEST(utils_enumerate_cuda, sequence_range) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // vecmem::cuda::device_memory_resource device_resource;
    // vecmem::cuda::host_memory_resource host_resource;

    vecmem::vector<dindex> reference({2, 3, 4, 5, 6, 7});

    const darray<dindex, 2> range = {static_cast<dindex>(2),
                                     static_cast<dindex>(7)};

    vecmem::data::vector_buffer<dindex> check_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(
            range[1] - range[0] + 1),
        0, managed_resource);

    /*
     vecmem::data::vector_buffer<dindex> check_buffer(
         static_cast<vecmem::data::vector_buffer<dindex>::size_type>(
             range[1] - range[0] + 1),
         0, device_resource);
     */
    copy.setup(check_buffer);

    sequence_range(range, check_buffer);

    vecmem::vector<dindex> check{&managed_resource};
    copy(check_buffer, check);

    ASSERT_EQ(check, reference);
}

// This tests the convenience enumeration function
TEST(utils_enumerate_cuda, enumerate) {

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}};

    for (auto [i, v] : enumerate(seq)) {
        ASSERT_EQ(i, v.ui);
    }
}

// This tests the restricted iterator
TEST(utils_enumerate_cuda, range) {
    size_t begin = 1;
    size_t end = 4;

    dvector<int> seq = {0, 1, 2, 3, 4, 5};

    size_t i = 1;
    for (const auto &v :
         iterator_range(seq, std::array<size_t, 2>{begin, end})) {
        ASSERT_NE(v, 5);
        ASSERT_EQ(v, seq[i++]);
    }
}
