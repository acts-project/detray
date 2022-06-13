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

void get_sum(tuple_vector_container_data_type& container_data,
             vecmem::data::vector_view<double>& sum_data);

}  // namespace detray