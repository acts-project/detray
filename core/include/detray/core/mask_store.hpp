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

/** Tuple based container with vector elements.
 *
 * @tparam tuple_t is the type of tuple
 * @tparam vector_t is the type of vector
 * @tparam id_t is the type of indexing integer
 * @tparam mask_types are the types of masks
 */
template <template <typename...> class tuple_t = dtuple,
          template <typename> class vector_t = dvector,
          typename id_t = unsigned int, typename... mask_types>
class mask_store final
    : public tuple_vector_container<tuple_t, vector_t, id_t, mask_types...> {

    public:
    // Convenient type declarations
    using base_type =
        tuple_vector_container<tuple_t, vector_t, id_t, mask_types...>;
    using base_type::base_type;
    using id_type = typename base_type::id_type;
    using link_type = typename base_type::link_type;
    using range_type = typename base_type::range_type;
    using mask_tuple = typename base_type::container_type;

    /**
     * Constructor with a vecmem memory resource. (host-side only)
     *
     * @param resource is the vecmem memory resource
     */
    DETRAY_HOST
    mask_store(vecmem::memory_resource &resource) : base_type(resource) {}

    /**
     * Constructor with a data container.
     * Mainly used in the device side
     *
     * @tparam is the type of input data container
     *
     * @param container_data is the data container
     */
    template <typename container_data_t>
    DETRAY_DEVICE mask_store(container_data_t &container_data)
        : base_type(container_data) {}

    template <id_type mask_id, typename... bounds_type>
    DETRAY_HOST auto &add_mask(bounds_type &&... mask_bounds) noexcept(false) {
        return base_type::template add_value<mask_id, bounds_type...>(
            std::forward<bounds_type>(mask_bounds)...);
    }

    template <std::size_t current_id = 0, typename mask_type>
    DETRAY_HOST inline void add_masks(vector_t<mask_type> &masks) noexcept(
        false) {
        return base_type::template add_vector<current_id, mask_type>(masks);
    }

    template <std::size_t current_id = 0, typename mask_type>
    DETRAY_HOST inline void add_masks(vector_t<mask_type> &&masks) noexcept(
        false) {
        return base_type::template add_vector<current_id, mask_type>(masks);
    }

    template <std::size_t current_id = 0>
    DETRAY_HOST inline void append_masks(mask_store &&other) {
        return base_type::template add_container<current_id>(other);
    }
};

// Alias of tuple_vector_container_data
template <typename mask_store_t>
using mask_store_data = tuple_vector_container_data<mask_store_t>;

}  // namespace detray
