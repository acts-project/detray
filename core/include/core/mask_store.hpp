/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "definitions/basic_types.hpp"
#include "definitions/detray_qualifiers.hpp"
#include "utils/enumerate.hpp"
#include "utils/expand.hpp"
#include "utils/indexing.hpp"

namespace detray {

/** A mask store that provides the correct mask containers to client classes. */
template <template <typename...> class vector_type = dvector,
          typename... mask_types>
class mask_store {

    public:
    using mask_tuple = __tuple::tuple<vector_type<mask_types>...>;

    /** Default constructor **/
    mask_store() = default;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    mask_store(vecmem::memory_resource &resource)
        : _mask_tuple(vector_type<mask_types>{&resource}...) {}

    /** Constructor with mask_store_data **/
#if defined(__CUDACC__)
    template <typename mask_store_data_t>
    DETRAY_DEVICE mask_store(mask_store_data_t &store_data)
        : _mask_tuple(
              store_data.device(typename gens<sizeof...(mask_types)>::type())) {
    }
#endif

    /** Size : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return the size of the vector containing the masks of the required type
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE size_t size() const {
        return __tuple::get<mask_id>(_mask_tuple).size();
    }

    /** Empty : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return whether the vector containing the masks of the required type
     * is empty
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE bool empty() const {
        return __tuple::get<mask_id>(_mask_tuple).empty();
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
    template <unsigned int current_id = 0>
    DETRAY_HOST_DEVICE const auto &mask(const dindex mask_id,
                                        const dindex mask_index) const {
        if (current_id == mask_id) {
            return group<current_id>()[mask_index];
        }
        // Next mask type
        if constexpr (current_id < __tuple::tuple_size<mask_tuple>::value - 1) {
            return mask<current_id + 1>(mask_id, mask_index);
        }
    }

    /** Retrieve a vector of masks of a certain type (mask group)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE auto &group() {
        return __tuple::get<mask_id>(_mask_tuple);
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE constexpr const auto &group() const {
        return __tuple::get<mask_id>(_mask_tuple);
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
    const auto &masks() { return _mask_tuple; }

    /** Obtain the vecmem data of mask store
     *
     * @return tuple type of vecmem::data::vector_view objects
     */
    template <unsigned int... S>
    DETRAY_HOST __tuple::tuple<vecmem::data::vector_view<mask_types>...> data(
        seq<S...>) {
        return std::make_tuple(vecmem::data::vector_view<mask_types>(
            vecmem::get_data(__tuple::get<S>(_mask_tuple)))...);
    }

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
        auto &mask_group = __tuple::get<mask_id>(_mask_tuple);
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
        auto &mask_group = __tuple::get<current_id>(_mask_tuple);

        if constexpr (std::is_same_v<decltype(masks), decltype(mask_group)>) {
            // Reserve memory and copy new masks
            mask_group.reserve(mask_group.size() + masks.size());
            mask_group.insert(mask_group.end(), masks.begin(), masks.end());
        }

        // Next mask type
        if constexpr (current_id < std::tuple_size_v<mask_tuple> - 1) {
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
        auto &mask_group = __tuple::get<current_id>(_mask_tuple);

        if constexpr (std::is_same_v<decltype(masks), decltype(mask_group)>) {
            // Reserve memory and copy new masks
            mask_group.reserve(mask_group.size() + masks.size());
            mask_group.insert(mask_group.end(),
                              std::make_move_iterator(masks.begin()),
                              std::make_move_iterator(masks.end()));
        }

        // Next mask type
        if constexpr (current_id < std::tuple_size_v<mask_tuple> - 1) {
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
        auto &mask_group = __tuple::get<current_id>(other);
        add_masks(mask_group);

        // Next mask type
        if constexpr (current_id < std::tuple_size_v<mask_tuple> - 1) {
            return append_masks<current_id + 1>(other);
        }
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
template <typename... mask_types>
struct mask_store_data {

    /** Constructor from mask store
     *
     * @param store is the mask store from host
     **/
    template <template <typename...> class vector_type>
    mask_store_data(mask_store<vector_type, mask_types...> &store)
        : _data(store.data(typename gens<sizeof...(mask_types)>::type())) {}

    /** Size : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return the size of the vector containing the masks of the required type
     */
    template <unsigned int mask_id>
    size_t size() const {
        return std::get<mask_id>(_data).size();
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    const auto &group() const {
        return std::get<mask_id>(_data);
    }

    /** Obtain the vecmem device vector of mask store
     *
     * @return tuple type of vecmem::data::vector_view objects
     */
    template <unsigned int... S>
    DETRAY_DEVICE __tuple::tuple<vecmem::device_vector<mask_types>...> device(
        seq<S...>) {
        return __tuple::make_tuple(
            vecmem::device_vector<mask_types>(std::get<S>(_data))...);
    }

    /** tuple of vecmem data **/
    std::tuple<vecmem::data::vector_view<mask_types>...> _data;
};

/** Get transform_store_data
 **/
template <template <typename...> class vector_type, typename... mask_types>
inline mask_store_data<mask_types...> get_data(
    mask_store<vector_type, mask_types...> &store) {
    return mask_store_data<mask_types...>(store);
}

}  // namespace detray
