/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray Test include(s)
#include "tests/common/test_defs.hpp"

// Detray Core include(s)
#include "detray/core/detail/tuple_array_container.hpp"
#include "detray/core/detail/tuple_vector_container.hpp"

// Vecmem include(s)
#include "vecmem/containers/device_vector.hpp"

namespace detray {

using tuple_vector_container_type =
    tuple_vector_container<thrust::tuple, vecmem::vector, std::size_t, int,
                           float, double>;

using tuple_vector_container_data_type =
    tuple_vector_container_data<tuple_vector_container_type>;

template <template <typename> class vector_t>
struct int_type {
    using object_type = vector_t<int>;
    using data_type = vecmem::data::vector_view<int>;
    using view_type = vecmem::data::vector_view<int>;
    static constexpr std::size_t N = 1;
};

template <template <typename> class vector_t>
struct float_type {
    using object_type = vector_t<float>;
    using data_type = vecmem::data::vector_view<float>;
    using view_type = vecmem::data::vector_view<float>;
    static constexpr std::size_t N = 2;
};

using tuple_array_container_type =
    tuple_array_container<thrust::tuple, std::array, std::size_t,
                          int_type<vecmem::vector>, float_type<vecmem::vector>>;

using tuple_array_container_view_type =
    tuple_array_container_view<tuple_array_container_type>;

void get_sum(tuple_vector_container_data_type& container_data,
             vecmem::data::vector_view<double>& sum_data);

void get_sum(tuple_array_container_view_type container_data,
             vecmem::data::vector_view<double>& sum_data);

}  // namespace detray