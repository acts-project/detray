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

#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/core/detail/tuple_container.hpp"

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

using tuple_t =
    detail::tuple_container<dtuple, vecmem::vector<int>, vecmem::vector<float>>;
using tuple_dev_t = detail::tuple_container<dtuple, vecmem::device_vector<int>,
                                            vecmem::device_vector<float>>;

template <typename view_t>
__global__ void tuple_kernel(view_t view) {
    std::size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;

    tuple_dev_t store(view);

    if (globalIdx >= detail::get<0>(store).size()) {
        return;
    }

    printf("%d\n", detail::get<0>(store).at(globalIdx));
    printf("%f\n", detail::get<1>(store).at(globalIdx));
}

GTEST_TEST(detray_detail, tuple_container) {

    // Run the kernel by copying the data to device

    vecmem::host_memory_resource host_mr;
    tuple_t store(host_mr);
    detail::get<0>(store).push_back(3);
    detail::get<1>(store).push_back(4.f);

    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy cpy;

    auto store_buffer = detray::get_buffer(store, dev_mr, cpy);

    tuple_kernel<<<8, 1>>>(detray::get_data(store_buffer));

    // Run the same kernel with managed memory

    vecmem::cuda::managed_memory_resource mng_mr;
    tuple_t managed_store(mng_mr);
    detail::get<0>(managed_store).push_back(5);
    detail::get<1>(managed_store).push_back(6.f);

    tuple_kernel<<<8, 1>>>(detray::get_data(managed_store));
}

template <template <typename...> class vector_t = dvector>
struct test {

    using view_type = dmulti_view<dvector_view<int>, dvector_view<float>>;
    using buffer_type =
        dmulti_buffer<dvector_buffer<int>, dvector_buffer<float>>;

    test() = default;
    template <typename view_t>
    test(view_t v)
        : first(detail::get<0>(v.m_view)), second(detail::get<1>(v.m_view)) {}

    view_type get_data() {
        return {vecmem::get_data(first), vecmem::get_data(second)};
    }

    vector_t<int> first;
    vector_t<float> second;
};

using multi_store_t = multi_store<std::size_t, empty_context, dtuple,
                                  vecmem::vector<float>, test<vecmem::vector>>;
using multi_store_dev_t =
    multi_store<std::size_t, empty_context, dtuple,
                vecmem::device_vector<float>, test<vecmem::device_vector>>;

template <typename view_t>
__global__ void multi_kernel(view_t view) {
    std::size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;

    multi_store_dev_t store(view);

    if (globalIdx >= store.get<0>().size()) {
        return;
    }

    printf("%d\n", store.get<0>().at(globalIdx));
    // printf("%f\n", store.get<1>().at(globalIdx));
}

GTEST_TEST(detray_detail, multi_store) {

    // Run the kernel by copying the data to device

    vecmem::host_memory_resource host_mr;
    multi_store_t store(host_mr);
    store.get<0>().push_back(3);
    // store.get<1>().push_back(4.f);

    vecmem::cuda::device_memory_resource dev_mr;
    vecmem::cuda::copy cpy;

    test<> t{};
    t.first.push_back(42);
    t.second.push_back(43.f);
    test<>::view_type test_view = t.get_data();
    test<>::buffer_type test_buff = detray::get_buffer(t, dev_mr, cpy);
    test<vecmem::device_vector> t_dev(test_view);

    auto store_buffer = detray::get_buffer(store, dev_mr, cpy);

    // using bla = typename decltype(store_buffer)::blub;

    multi_kernel<<<8, 1>>>(detray::get_data(store_buffer));

    // Run the same kernel with managed memory

    vecmem::cuda::managed_memory_resource mng_mr;
    multi_store_t managed_store(mng_mr);
    managed_store.get<0>().push_back(5);
    // managed_store.get<1>().push_back(6.f);

    multi_kernel<<<8, 1>>>(detray::get_data(managed_store));
}