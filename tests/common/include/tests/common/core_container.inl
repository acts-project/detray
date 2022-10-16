/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detail/multi_type_store.hpp"
#include "detray/core/detail/new_tuple_container.hpp"
#include "detray/core/detail/tuple_array_container.hpp"
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
    output_type operator()(const container_t& coll, const int /*index*/) {
        return coll.size();
    }
};

TEST(container, tuple_container) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;
    using host_tuple_t =
        detail::tuple_container<std::tuple, vecmem::vector<float>,
                                vecmem::vector<std::size_t>,
                                vecmem::vector<int*>>;

    using device_tuple_t =
        detail::tuple_container<thrust::tuple, vecmem::device_vector<float>,
                                vecmem::device_vector<std::size_t>,
                                vecmem::device_vector<int*>>;

    vecmem::vector<float> vec1{};
    vecmem::vector<std::size_t> vec2{};
    vecmem::vector<int*> vec3{};

    host_tuple_t container(resource, vec1, vec2, vec3);

    static_assert(
        std::is_same_v<
            typename host_tuple_t::view_type,
            dmulti_view<dvector_view<float>, dvector_view<std::size_t>,
                        dvector_view<int*>>>,
        "View type incorrectly assembled");

    static_assert(std::is_same_v<typename host_tuple_t::const_view_type,
                                 dmulti_view<dvector_view<const float>,
                                             dvector_view<const std::size_t>,
                                             dvector_view<int* const>>>,
                  "Const view type incorrectly assembled");

    typename host_tuple_t::view_type view = get_data(container);
    device_tuple_t dev_container(view);

    // Base container function check
    EXPECT_EQ(container.size(), 3);
    EXPECT_EQ(dev_container.size(), 3);

    EXPECT_TRUE(container.get<0>().empty());
    EXPECT_TRUE(container.get<1>().empty());
    EXPECT_TRUE(detail::get<2>(container).empty());
}

TEST(container, vector_multi_type_store) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;

    // Create tuple vector container
    multi_type_store<std::size_t, empty_context, std::tuple, vecmem::vector,
                     std::size_t, float, double>
        vector_store(resource);

    // Base container function check
    EXPECT_EQ(vector_store.n_collections(), 3);
    EXPECT_EQ(vector_store.empty<0>(), true);
    EXPECT_EQ(vector_store.empty<1>(), true);
    EXPECT_EQ(vector_store.empty<2>(), true);
    EXPECT_EQ(vector_store.size<0>(), 0UL);
    EXPECT_EQ(vector_store.size<1>(), 0UL);
    EXPECT_EQ(vector_store.size<2>(), 0UL);

    // Add elements to the container
    vector_store.push_back<0>(1UL);
    vector_store.emplace_back<0>(empty_context{}, 2UL);
    vector_store.push_back<1>(3.1f);
    vector_store.emplace_back<1>(empty_context{}, 4.5f);
    vector_store.push_back<2>(5.5);
    vector_store.emplace_back<2>(empty_context{}, 6.);

    EXPECT_EQ(vector_store.empty<0>(), false);
    EXPECT_EQ(vector_store.empty<1>(), false);
    EXPECT_EQ(vector_store.empty<2>(), false);
    EXPECT_EQ(vector_store.size<0>(), 2UL);
    EXPECT_EQ(vector_store.size<1>(), 2UL);
    EXPECT_EQ(vector_store.size<2>(), 2UL);

    vecmem::vector<std::size_t> int_vec{3UL, 4UL, 5UL};
    vector_store.insert(int_vec);

    vecmem::vector<float> float_vec{12.1f};
    vector_store.insert(float_vec);

    vector_store.insert(vecmem::vector<double>{10.5, 7.6});

    // int collectiont
    EXPECT_EQ(vector_store.size<0>(), 5UL);
    EXPECT_EQ(vector_store.get<0>()[0], 1UL);
    EXPECT_EQ(vector_store.get<0>()[1], 2UL);
    EXPECT_EQ(vector_store.get<0>()[2], 3UL);
    EXPECT_EQ(vector_store.get<0>()[3], 4UL);
    EXPECT_EQ(vector_store.get<0>()[4], 5UL);

    // float collectiont
    EXPECT_EQ(vector_store.size<1>(), 3UL);
    EXPECT_FLOAT_EQ(vector_store.get<1>()[0], 3.1f);
    EXPECT_FLOAT_EQ(vector_store.get<1>()[1], 4.5f);
    EXPECT_FLOAT_EQ(vector_store.get<1>()[2], 12.1f);

    // double collectiont
    EXPECT_EQ(vector_store.size<2>(), 4UL);
    EXPECT_FLOAT_EQ(vector_store.get<2>()[0], 5.5);
    EXPECT_FLOAT_EQ(vector_store.get<2>()[1], 6.);
    EXPECT_FLOAT_EQ(vector_store.get<2>()[2], 10.5);
    EXPECT_FLOAT_EQ(vector_store.get<2>()[3], 7.6);

    // call functor
    EXPECT_EQ(vector_store.call<test_func>(std::make_pair(0, 0)), 5UL);
    EXPECT_EQ(vector_store.call<test_func>(std::make_pair(1, 0)), 3UL);
    EXPECT_EQ(vector_store.call<test_func>(std::make_pair(2, 0)), 4UL);
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
