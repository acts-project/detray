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
          template <typename...> class vector_t = dvector,
          typename ID = unsigned int, typename... materials_t>
struct material_store {

    public:
    template <typename... Args>
    using tuple_type = tuple_t<Args...>;

    /**
     * material_tuple does not follow the tuple_type.
     * vtuple has different types based on the file location 1) std::tuple in
     * cpp/hpp; 2) thrust::tuple in cu
     */
    using material_tuple = vtuple::tuple<vector_t<materials_t>...>;

    /** data type for material_store_data **/
    using material_tuple_data =
        tuple_type<vecmem::data::vector_view<mask_types>...>;

    /** Default constructor **/
    material_store() = delete;

    /**
     * Constructor with vecmem memory resource
     *
     * @param resource is the vecmem memory resource
     */
    DETRAY_HOST
    material_store(vecmem::memory_resource &resource)
        : _material_tuple(vector_t<materials_t>{&resource}...) {}

    /**
     * Constructor with mask_store_data
     *
     * @param store_data is the material store data
     */
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