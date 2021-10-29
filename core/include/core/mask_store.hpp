/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "definitions/detray_get.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {

/** A mask store that provides the correct mask containers to client classes. */
template <template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          typename... mask_types>
class mask_store {

    public:
    // using mask_tuple = tuple_type<vector_type<mask_types>...>;
    using mask_tuple = dtuple<vector_type<mask_types>...>;

    /** Size : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return the size of the vector containing the masks of the required type
     */
    template <unsigned int mask_id>
    const size_t size() const {
        return detail::get<mask_id>(_mask_tuple).size();
    }

    /** Empty : Contextual STL like API
     *
     * @tparam mask_id the index for the mask_type
     * @return whether the vector containing the masks of the required type
     * is empty
     */
    template <unsigned int mask_id>
    bool empty() const {
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
    template <unsigned int current_id = 0>
    const auto &mask(const unsigned int mask_id,
                     const dindex mask_index) const {
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
    constexpr auto &group() {
        return detail::get<mask_id>(_mask_tuple);
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <unsigned int mask_id>
    constexpr const auto &group() const {
        return detail::get<mask_id>(_mask_tuple);
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
    inline void add_masks(vector_type<mask_type> &masks) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = detail::get<current_id>(_mask_tuple);

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
    inline void add_masks(vector_type<mask_type> &&masks) noexcept(false) {
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
    inline void append_masks(mask_store &&other) {
        // Add masks to current group
        auto &mask_group = detail::get<current_id>(other);
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

}  // namespace detray
