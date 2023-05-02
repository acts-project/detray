/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray Core include(s)
#include "detray/core/detail/multi_store.hpp"

// Vecmem include(s)
#include "vecmem/containers/device_vector.hpp"

// Thrust include(s)
#include <thrust/tuple.h>

namespace detray {

using host_store_type =
    regular_multi_store<int, empty_context, dtuple, vecmem::vector, std::size_t,
                        float, double>;

using device_store_type =
    regular_multi_store<int, empty_context, dtuple, vecmem::device_vector,
                        std::size_t, float, double>;

void get_sum(typename host_store_type::view_type store_view,
             vecmem::data::vector_view<double> sum_data);

}  // namespace detray
