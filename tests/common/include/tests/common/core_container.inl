/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detail/tuple_array_container.hpp"
#include "detray/core/detail/tuple_vector_container.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"

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

struct test_func {
    using output_type = std::size_t;

    template <typename container_t>
    output_type operator()(const container_t& gr, const int /*index*/) {
        return gr.size();
    }
};

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

    vecmem::vector<float> float_vec{12.1};
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
    EXPECT_EQ(container.size<1>(), 3);
    EXPECT_FLOAT_EQ(container.group<1>()[0], 3.1);
    EXPECT_FLOAT_EQ(container.group<1>()[1], 4.5);
    EXPECT_FLOAT_EQ(container.group<1>()[2], 12.1);

    // double group
    EXPECT_EQ(container.size<2>(), 4);
    EXPECT_FLOAT_EQ(container.group<2>()[0], 5.5);
    EXPECT_FLOAT_EQ(container.group<2>()[1], 6.);
    EXPECT_FLOAT_EQ(container.group<2>()[2], 10.5);
    EXPECT_FLOAT_EQ(container.group<2>()[3], 7.6);

    // unrolling test
    EXPECT_EQ(container.call<test_func>(std::make_pair(0, 0)), 5);
    EXPECT_EQ(container.call<test_func>(std::make_pair(1, 0)), 3);
    EXPECT_EQ(container.call<test_func>(std::make_pair(2, 0)), 4);
}

using grid2r =
    grid2<replace_populator, axis::regular, axis::regular, serializer2>;

using grid2a =
    grid2<attach_populator, axis::regular, axis::regular, serializer2>;

TEST(container, tuple_array_container) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;

    // Create tuple array container
    tuple_array_container<std::tuple, std::array, std::size_t,
                          std::index_sequence<1, 2>, grid2r, grid2a>
        container(resource);

    EXPECT_EQ(container.size(), 2);
    EXPECT_EQ(container.empty<0>(), false);
    EXPECT_EQ(container.empty<1>(), false);
    EXPECT_EQ(container.to_id<>(0), 0);
    EXPECT_EQ(container.to_id<1>(1), 1);
    EXPECT_EQ(container.to_id<1>(0), 2);

    // Populate the elements
    auto& grid2_replace_array = container.group<0>();
    grid2r::axis_p0_type xaxisr{1, 0., 1., resource};
    grid2r::axis_p1_type yaxisr{1, 0., 1., resource};
    grid2_replace_array[0] =
        grid2r(std::move(xaxisr), std::move(yaxisr), resource);
    grid2_replace_array[0].bin(0, 0) = 3;

    auto& grid2_attach_array = container.group<1>();
    grid2a::axis_p0_type xaxisa{2, 0., 2., resource};
    grid2a::axis_p1_type yaxisa{1, 0., 1., resource};

    grid2_attach_array[0] =
        grid2a(std::move(xaxisa), std::move(yaxisa), resource);
    grid2_attach_array[0].bin(0, 0).push_back(0);
    grid2_attach_array[0].bin(0, 0).push_back(1);
    grid2_attach_array[0].bin(1, 0).push_back(2);
    grid2_attach_array[0].bin(1, 0).push_back(3);

    grid2_attach_array[1] =
        grid2a(std::move(xaxisa), std::move(yaxisa), resource);
    grid2_attach_array[1].bin(0, 0).push_back(4);
    grid2_attach_array[1].bin(1, 0).push_back(5);

    // Check the element
    EXPECT_EQ(container.size<0>(), 1);
    EXPECT_EQ(container.group<0>()[0].bin(0, 0), 3);

    EXPECT_EQ(container.size<1>(), 2);
    EXPECT_EQ(container.group<1>()[0].bin(0, 0),
              vecmem::vector<dindex>({0, 1}));
    EXPECT_EQ(container.group<1>()[0].bin(1, 0),
              vecmem::vector<dindex>({2, 3}));
    EXPECT_EQ(container.group<1>()[1].bin(0, 0), vecmem::vector<dindex>({4}));
    EXPECT_EQ(container.group<1>()[1].bin(1, 0), vecmem::vector<dindex>({5}));
}
