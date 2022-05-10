/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detail/tuple_array_container.hpp"
#include "detray/core/detail/tuple_vector_container.hpp"

// Vecmem include(s)
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <tuple>
#include <vector>

using namespace detray;

TEST(container, tuple_vector_container) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;

    // Create tuple vector container
    tuple_vector_container<std::tuple, vecmem::vector, std::size_t, int, float,
                           double>
        container(resource);

    // Base container function check
    EXPECT_EQ(container.size(), 3);
    EXPECT_EQ(container.empty<0>(), true);
    EXPECT_EQ(container.empty<1>(), true);
    EXPECT_EQ(container.empty<2>(), true);

    EXPECT_EQ(container.to_id<>(0), 0);
    EXPECT_EQ(container.to_id<1>(2), 2);
    EXPECT_EQ(container.to_id<1>(0), 3);

    // Add elements to the container
    container.add_value<0>(1);
    container.add_value<0>(2);
    container.add_value<1>(3.1);
    container.add_value<1>(4.5);
    container.add_value<2>(5.5);
    container.add_value<2>(6.);

    vecmem::vector<int> int_vec{3, 4, 5};
    container.add_vector(int_vec);

    vecmem::vector<float> float_vec{12.1, 5.6};
    container.add_vector(float_vec);

    container.add_vector(vecmem::vector<double>{10.5, 7.6});

    EXPECT_EQ(container.empty<0>(), false);
    EXPECT_EQ(container.empty<1>(), false);
    EXPECT_EQ(container.empty<2>(), false);

    // int group
    EXPECT_EQ(container.size<0>(), 5);
    EXPECT_EQ(container.group<0>()[0], 1);
    EXPECT_EQ(container.group<0>()[1], 2);
    EXPECT_EQ(container.group<0>()[2], 3);
    EXPECT_EQ(container.group<0>()[3], 4);
    EXPECT_EQ(container.group<0>()[4], 5);

    // float group
    EXPECT_EQ(container.size<1>(), 4);
    EXPECT_FLOAT_EQ(container.group<1>()[0], 3.1);
    EXPECT_FLOAT_EQ(container.group<1>()[1], 4.5);
    EXPECT_FLOAT_EQ(container.group<1>()[2], 12.1);
    EXPECT_FLOAT_EQ(container.group<1>()[3], 5.6);

    // double group
    EXPECT_EQ(container.size<2>(), 4);
    EXPECT_FLOAT_EQ(container.group<2>()[0], 5.5);
    EXPECT_FLOAT_EQ(container.group<2>()[1], 6.);
    EXPECT_FLOAT_EQ(container.group<2>()[2], 10.5);
    EXPECT_FLOAT_EQ(container.group<2>()[3], 7.6);
}

struct int_type {
    using object_type = vecmem::vector<int>;
    using data_type = vecmem::data::vector_view<int>;
    using view_type = vecmem::data::vector_view<int>;
    static constexpr std::size_t N = 1;
};

struct float_type {
    using object_type = vecmem::vector<float>;
    using data_type = vecmem::data::vector_view<float>;
    using view_type = vecmem::data::vector_view<float>;
    static constexpr std::size_t N = 2;
};

TEST(container, tuple_array_container) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;

    // Create tuple array container
    tuple_array_container<std::tuple, std::array, std::size_t, int_type,
                          float_type>
        container(resource);

    // Base container function check
    EXPECT_EQ(container.size(), 2);
    EXPECT_EQ(container.empty<0>(), false);
    EXPECT_EQ(container.empty<1>(), false);
    EXPECT_EQ(container.to_id<>(0), 0);
    EXPECT_EQ(container.to_id<1>(1), 1);
    EXPECT_EQ(container.to_id<1>(0), 2);

    // Populate the elements
    auto& int_group = container.group<0>();
    int_group = std::array<vecmem::vector<int>, 1>{
        vecmem::vector<int>{{3, 4, 5}, &resource}};

    auto& float_group = container.group<1>();
    float_group = std::array<vecmem::vector<float>, 2>{
        vecmem::vector<float>{{1.1, 4.2}, &resource},
        vecmem::vector<float>{{5.1, 2.3, 1.7}, &resource}};

    // int group
    EXPECT_EQ(container.size<0>(), 1);
    EXPECT_EQ(container.group<0>()[0], vecmem::vector<int>({3, 4, 5}));

    // float group
    EXPECT_EQ(container.size<1>(), 2);
    EXPECT_EQ(container.group<1>()[0], vecmem::vector<float>({1.1, 4.2}));
    EXPECT_EQ(container.group<1>()[1], vecmem::vector<float>({5.1, 2.3, 1.7}));
}