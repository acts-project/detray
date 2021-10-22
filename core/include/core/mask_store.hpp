/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <thrust/iterator/zip_iterator.h>

#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "definitions/detray_qualifiers.hpp"
#include "utils/enumerate.hpp"
#include "utils/expand.hpp"
#include "utils/indexing.hpp"

namespace detray {

/** Forward declaration of mask_store_data **/
template <typename... mask_types>
struct mask_store_data;

/** A mask store that provides the correct mask containers to client classes. */
template <template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          typename... mask_types>
class mask_store {

    public:
    using mask_tuple = tuple_type<vector_type<mask_types>...>;

    /** Default constructor **/
    mask_store() = default;

    /** Constructor with vecmem memory resource **/
    mask_store(vecmem::memory_resource &resource)
        : _mask_tuple(vector_type<mask_types>{&resource}...) {}

    /** Constructor with mask_store_data **/
#if defined(__CUDACC__)
    DETRAY_DEVICE
    mask_store(mask_store_data<mask_types...> &store_data)
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
    const size_t size() const {
        return std::get<mask_id>(_mask_tuple).size();
    }

    /** Empty : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return whether the vector containing the masks of the required type
     * is empty
     */
    template <unsigned int mask_id>
    bool empty() const {
        return std::get<mask_id>(_mask_tuple).empty();
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
    const auto &mask(const dindex mask_id, const dindex mask_index) const {
        if (current_id == mask_id) {
            return group<current_id>()[mask_index];
        }
        // Next mask type
        if constexpr (current_id < std::tuple_size_v<mask_tuple> - 1) {
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
#if defined(__CUDACC__)
        return thrust::get<mask_id>(_mask_tuple);
#else
        return std::get<mask_id>(_mask_tuple);
#endif
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    DETRAY_HOST_DEVICE const auto &group() const {
#if defined(__CUDACC__)
        return thrust::get<mask_id>(_mask_tuple);
#else
        return std::get<mask_id>(_mask_tuple);
#endif
    }

    /** Access underlying container - const
     *
     * @return internal masks tuple type
     */
    const auto &masks() const { return _mask_tuple; }

    /** Access underlying container
     *
     * @return internal masks tuple type
     */
    const auto &masks() { return _mask_tuple; }

    /** Obtain the vecmem data of mask store
     *
     * @return tuple type of vecmem::data::vector_view objects
     */
    template <unsigned int... S>
    thrust::tuple<vecmem::data::vector_view<mask_types>...> data(seq<S...>) {
        return thrust::make_tuple(vecmem::data::vector_view<mask_types>(
            vecmem::get_data(std::get<S>(_mask_tuple)))...);
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
    auto &add_mask(bounds_type &&... mask_bounds) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_id>(_mask_tuple);
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
    template <unsigned int mask_id, typename mask_type>
    void add_masks(const vector_type<mask_type> &masks) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_id>(_mask_tuple);
        // Reserve memory and copy new masks
        mask_group.reserve(mask_group.size() + masks.size());
        mask_group.insert(mask_group.end(), masks.begin(), masks.end());
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
    template <unsigned int mask_id, typename mask_type>
    void add_masks(vector_type<mask_type> &&masks) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_id>(_mask_tuple);
        // Reserve memory and copy new masks
        mask_group.reserve(mask_group.size() + masks.size());
        mask_group.insert(mask_group.end(),
                          std::make_move_iterator(masks.begin()),
                          std::make_move_iterator(masks.end()));
    }

    /** Append a mask store to the current one
     *
     * @tparam current_index to start unrolling at
     *
     * @param other The other mask store, move semantics
     *
     * @note in general can throw an exception
     */
    template <unsigned int current_index = 0>
    inline void append_masks(mask_store &&other) {
        // Add masks to current group
        auto &mask_group = std::get<current_index>(other);
        add_masks<current_index>(mask_group);

        // Next mask type
        if constexpr (current_index < std::tuple_size_v<mask_tuple> - 1) {
            return append_masks<current_index + 1>(other);
        }
    }

    private:
    /** tuple of mask vectors (mask groups) */
    mask_tuple _mask_tuple;
};

/** A static inplementation of mask store data for device*/
template <typename... mask_types>
struct mask_store_data {

    /** Constructor from mask store
     *
     * @param store is the mask store from host
     **/
    template <template <typename...> class tuple_type,
              template <typename...> class vector_type>
    mask_store_data(mask_store<tuple_type, vector_type, mask_types...> &store)
        : _data(store.data(typename gens<sizeof...(mask_types)>::type())) {}

    /** Size : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return the size of the vector containing the masks of the required type
     */
    template <unsigned int mask_id>
    const size_t size() const {
        return thrust::get<mask_id>(_data).size();
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    const auto &group() const {
        return thrust::get<mask_id>(_data);
    }

    /** Obtain the vecmem device vector of mask store
     *
     * @return tuple type of vecmem::data::vector_view objects
     */
#if defined(__CUDACC__)
    template <unsigned int... S>
    DETRAY_DEVICE thrust::tuple<vecmem::device_vector<mask_types>...> device(
        seq<S...>) {
        return thrust::make_tuple(
            vecmem::device_vector<mask_types>(thrust::get<S>(_data))...);
    }
#endif

    /** tuple of vecmem data **/
    thrust::tuple<vecmem::data::vector_view<mask_types>...> _data;
};

/** Get transform_store_data
 **/
template <template <typename...> class tuple_type,
          template <typename...> class vector_type, typename... mask_types>
inline mask_store_data<mask_types...> get_data(
    mask_store<tuple_type, vector_type, mask_types...> &store) {
    return mask_store_data<mask_types...>(store);
}

}  // namespace detray
