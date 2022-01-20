/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>
#include <utility>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"
namespace detray {

/** A mask store that provides the correct mask containers to client classes. */
template <template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          typename... mask_types>
class mask_store {

    public:
    /**
     * mask_tuple is the only member which does not follow the tuple_type.
     * vtuple has different types based on the file location 1) std::tuple in
     * cpp/hpp; 2) thrust::tuple in cu
     */
    using mask_tuple = vtuple::tuple<vector_type<mask_types>...>;

    /** data type for mask_store_data **/
    using mask_tuple_data =
        tuple_type<vecmem::data::vector_view<mask_types>...>;

    // How to index data in this container
    /** Link type for a single mask: group type, index in group */
    using mask_link = std::array<dindex, 2>;
    /** Link type for a range of masks: group type, range in group */
    using mask_range = tuple_type<dindex, std::array<dindex, 2>>;

    /**
     * tuple_type for mask_tuple makes an illegal memory access error
     */
    // using mask_tuple = tuple_type<vector_type<mask_types>...>;

    /** Default constructor **/
    mask_store() = delete;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    mask_store(vecmem::memory_resource &resource)
        : _mask_tuple(vector_type<mask_types>{&resource}...) {}

    /** Constructor with mask_store_data **/
    template <
        typename mask_store_data_t,
        std::enable_if_t<
            !std::is_base_of_v<vecmem::memory_resource, mask_store_data_t> &&
                !std::is_same_v<mask_store, mask_store_data_t>,
            bool> = true>
    DETRAY_DEVICE mask_store(mask_store_data_t &store_data)
        : _mask_tuple(device(
              store_data, std::make_index_sequence<sizeof...(mask_types)>{})) {}

    /** Size : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return the size of the vector containing the masks of the required type
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE size_t size() const {
        return detail::get<mask_id>(_mask_tuple).size();
    }

    /** Size : Contextual STL like API
     * @return the number of mask types in the store
     */
    DETRAY_HOST_DEVICE constexpr size_t size() const {
        return detail::tuple_size<mask_tuple>::value;
    }

    /** Empty : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return whether the vector containing the masks of the required type
     * is empty
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE bool empty() const {
        return detail::get<mask_id>(_mask_tuple).empty();
    }

    /** Retrieve a single mask - const
     *
     * @tparam current_id the index for the mask_type
     *
     * @param mask_id the index for the mask_type
     * @param mask_index the index for the mask
     *
     * @return the required mask
     */
    /*template <unsigned int current_id = 0>
    DETRAY_HOST_DEVICE const auto &find_mask(const unsigned int mask_id,
                                        const dindex mask_index) const {
        if (current_id == mask_id) {
            return group<current_id>()[mask_index];
        }
        // Next mask type
        if constexpr (current_id < detail::tuple_size<mask_tuple>::value - 1) {
            return find_mask<current_id + 1>(mask_id, mask_index);
        }
    }*/

    /** Retrieve a single mask - const
     *
     * @tparam current_id the index for the mask_type
     *
     * @param mask_id the index for the mask_type
     * @param mask_index the index for the mask
     *
     * @return the required mask
     */
    /*DETRAY_HOST_DEVICE const auto &find_mask(mask_link m_link) const {
        return find_mask(detail::get<0>(m_link), detail::get<1>(m_link));
    }*/

    /** Retrieve a vector of masks of a certain type (mask group)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE constexpr auto &group() {
        return detail::get<mask_id>(_mask_tuple);
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE constexpr const auto &group() const {
        return detail::get<mask_id>(_mask_tuple);
    }

    /** Access underlying container - const
     *
     * @return internal masks tuple type
     */
    DETRAY_HOST_DEVICE
    const auto &masks() const { return _mask_tuple; }

    /** Access underlying container
     *
     * @return internal masks tuple type
     */
    DETRAY_HOST_DEVICE
    auto &masks() { return _mask_tuple; }

    /** Add a new mask in place
     *
     * @tparam mask_id index for this mask type in masks container
     * @tparam bounds_type type of the masks bounds
     *
     * @param mask_bounds list of mask bounds for construction
     *
     * @note in general can throw an exception
     */
    template <unsigned int mask_id, typename... bounds_type>
    DETRAY_HOST auto &add_mask(bounds_type &&... mask_bounds) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = detail::get<mask_id>(_mask_tuple);
        // Construct new mask in place
        return mask_group.emplace_back(
            std::forward<bounds_type>(mask_bounds)...);
    }

    /** Add a new bunch of masks
     *
     * @tparam mask_id index for this mask type in masks container
     * @tparam mask_type mask type to be updated
     *
     * @param masks Vector of masks to be added
     *
     * @note in general can throw an exception
     */
    template <unsigned int current_id = 0, typename mask_type>
    DETRAY_HOST inline void add_masks(vector_type<mask_type> &masks) noexcept(
        false) {
        // Get the mask group that will be updated
        auto &mask_group = detail::get<current_id>(_mask_tuple);

        if constexpr (std::is_same_v<decltype(masks), decltype(mask_group)>) {
            // Reserve memory and copy new masks
            mask_group.reserve(mask_group.size() + masks.size());
            mask_group.insert(mask_group.end(), masks.begin(), masks.end());
        }

        // Next mask type
        if constexpr (current_id < detail::tuple_size<mask_tuple>::value - 1) {
            return add_masks<current_id + 1>(masks);
        }
    }

    /** Add a new bunch of masks - move semantics
     *
     * @tparam mask_id index for this mask type in masks container
     * @tparam mask_type mask type to be updated
     *
     * @param masks Vector of masks to be added
     *
     * @note in general can throw an exception
     */
    template <unsigned int current_id = 0, typename mask_type>
    DETRAY_HOST inline void add_masks(vector_type<mask_type> &&masks) noexcept(
        false) {
        // Get the mask group that will be updated
        auto &mask_group = detail::get<current_id>(_mask_tuple);

        if constexpr (std::is_same_v<decltype(masks), decltype(mask_group)>) {
            // Reserve memory and copy new masks
            mask_group.reserve(mask_group.size() + masks.size());
            mask_group.insert(mask_group.end(),
                              std::make_move_iterator(masks.begin()),
                              std::make_move_iterator(masks.end()));
        }

        // Next mask type
        if constexpr (current_id < detail::tuple_size<mask_tuple>::value - 1) {
            return add_masks<current_id + 1>(masks);
        }
    }

    /** Append a mask store to the current one
     *
     * @tparam current_index to start unrolling at (if the mask id is known,
     *         unrolling can be started there)
     *
     * @param other The other mask store, move semantics
     *
     * @note in general can throw an exception
     */
    template <unsigned int current_id = 0>
    DETRAY_HOST inline void append_masks(mask_store &&other) {
        // Add masks to current group
        auto &mask_group = detail::get<current_id>(other);
        add_masks(mask_group);

        // Next mask type
        if constexpr (current_id < detail::tuple_size<mask_tuple>::value - 1) {
            return append_masks<current_id + 1>(other);
        }
    }

    /**
     * Get vecmem::data::vector_view objects
     */
    template <std::size_t... ints>
    DETRAY_HOST mask_tuple_data data(std::index_sequence<ints...> /*seq*/) {
        return detail::make_tuple<tuple_type>(
            vecmem::data::vector_view<mask_types>(
                vecmem::get_data(detail::get<ints>(_mask_tuple)))...);
    }

    /**
     * Get vecmem::device_vector objects
     */
    template <typename mask_store_data_t, std::size_t... ints>
    DETRAY_DEVICE mask_tuple device(mask_store_data_t &data,
                                    std::index_sequence<ints...> /*seq*/) {
        return vtuple::make_tuple(
            vector_type<mask_types>(detail::get<ints>(data._data))...);
    }

    private:
    /** tuple of mask vectors (mask groups) */
    mask_tuple _mask_tuple;
};

/** A static inplementation of mask store data for device
 *
 * The tuple type for mask store data is fixed to std::tuple since there was an
 * issue that _data gets corrupted when passed to .cu file
 */
template <typename mask_store_t>
struct mask_store_data {

    using mask_tuple_data_t = typename mask_store_t::mask_tuple_data;

    /** Constructor from mask store
     *
     * @param store is the mask store from host
     **/
    template <std::size_t... ints>
    DETRAY_HOST mask_store_data(mask_store_t &store,
                                std::index_sequence<ints...> seq)
        : _data(store.data(seq)) {}

    /** tuple of vecmem data
     */
    mask_tuple_data_t _data;
};

/** Get mask_store_data
 **/
template <template <typename...> class tuple_type,
          template <typename...> class vector_type, typename... mask_types>
inline mask_store_data<mask_store<tuple_type, vector_type, mask_types...>>
get_data(mask_store<tuple_type, vector_type, mask_types...> &store) {
    return mask_store_data<mask_store<tuple_type, vector_type, mask_types...>>(
        store, std::make_index_sequence<sizeof...(mask_types)>{});
}

}  // namespace detray
