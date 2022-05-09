/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Detray include(s)
#include "detray/core/detail/tuple_vector_container.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

// Mask store type declarations
template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_t,
          typename... mask_types>
using mask_store =
    tuple_vector_container<tuple_t, vector_t, id_t, mask_types...>;

template <typename mask_store_t>
using mask_store_data = tuple_vector_container_data<mask_store_t>;

}  // namespace detray
