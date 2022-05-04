/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/indexing.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

/** A material store that provides the correct material containers to client
 * classes. */
template <template <typename...> class tuple_t = dtuple,
          template <typename, std::size_t> class array_t = darray,
          typename ID = unsigned int, typename... Ts>
struct tuple_array_container {

    public:
    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;

    template <typename... Args>
    using tuple_type = tuple_t<Args...>;

    using container_tuple = tuple_type<array_type<Ts::object_type, Ts::N>...>;

    material_store() = delete;

    DETRAY_HOST
    material_store(vecmem::memory_resource &resource)
        : _material_tuple(vector_t<materials_t>{&resource}...) {}

    template <typename material_store_data_t,
              std::enable_if_t<
                  !std::is_base_of_v<vecmem::memory_resource,
                                     material_store_data_t> &&
                      !std::is_same_v<material_store, material_store_data_t>,
                  bool> = true>
    DETRAY_DEVICE mask_store(material_store_data_t &store_data)
        : _material_tuple(device(
              store_data, std::make_index_sequence<sizeof...(materials_t)>{})) {
    }
};

}  // namespace detray