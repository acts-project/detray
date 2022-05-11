/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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

// Material store type declarations
template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_t,
          typename... material_types>
using material_store =
    tuple_vector_container<tuple_t, vector_t, id_t, material_types...>;

template <typename material_store_t>
using material_store_data = tuple_vector_container_data<material_store_t>;

}  // namespace detray