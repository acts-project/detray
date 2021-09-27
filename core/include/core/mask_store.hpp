/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {

/** A mask store that provides the correct mask containers to client classes. */
template <template <typename...> class tuple_type = dtuple,
          template <typename> class vector_type = dvector,
          typename ...mask_types>
class mask_store {

    public:

    using mask_tuple = tuple_type<vector_type<mask_types>...>;

    /** Size : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    constexpr auto n_types() const {
        return std::tuple_size_v<mask_tuple>;
    }

    /** Size : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    template<unsigned int mask_index>
    const size_t size() const {
        return std::get<mask_index>(_mask_tuple).size();
    }

    /** Empty : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    template<unsigned int mask_index>
    bool empty() const {
        return std::get<mask_index>(_mask_tuple).empty();
    }

    /** Retrieve a vector of masks of a certain type (mask group)
     *
     * @tparam mask_index index of requested mask type in masks container
     *
     * @return vector of masks of a given type.
     */
    template<unsigned int mask_index>
    auto& group() {
        return std::get<mask_index>(_mask_tuple);
    }

    /** Retrieve a vector of masks of a certain type (mask group) - const
     *
     * @tparam mask_index index of requested mask type in masks container
     *
     * @return vector of masks of a given type.
     */
    template<unsigned int mask_index>
    const auto& group() const {
        return std::get<mask_index>(_mask_tuple);
    }

    /** Access underlying container - const
     *
     * @return internal masks tuple type
     */
    const auto& masks() const { return _mask_tuple; }

    /** Access underlying container
     *
     * @return internal masks tuple type
     */
    const auto& masks() { return _mask_tuple; }

    /** Add a new mask in place
     *
     * @tparam mask_index index for this mask type in masks container
     * @tparam bounds_type type of the masks bounds
     *
     * @param mask_bounds list of mask bounds for construction
     *
     * @note in general can throw an exception
     */
    template<unsigned int mask_index,
             typename mask_type>
    void add_mask(mask_type &&mask) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_index>(_mask_tuple);
        // Construct new mask in place
        mask_group.push_back(std::forward<mask_type>(mask));
    }

    /** Add a new mask in place
     *
     * @tparam mask_index index for this mask type in masks container
     * @tparam bounds_type type of the masks bounds
     *
     * @param mask_bounds list of mask bounds for construction
     *
     * @note in general can throw an exception
     */
    template<unsigned int mask_index,
             typename... bounds_type>
    void add_mask(bounds_type &&...mask_bounds) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_index>(_mask_tuple);
        // Construct new mask in place
        mask_group.emplace_back(std::forward<bounds_type>(mask_bounds)...);
    }

    /** Add a new bunch of masks
     *
     * @tparam mask_index index for this mask type in masks container
     * @tparam mask_type mask type to be updated
     *
     * @param masks Vector of masks to be added
     *
     * @note in general can throw an exception
     */
    template<unsigned int mask_index,
             typename mask_type>
    void add_masks(const vector_type<mask_type> &masks) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_index>(_mask_tuple);
        // Reserve memory and copy new masks
        mask_group.reserve(mask_group.size() + masks.size());
        mask_group.insert(mask_group.end(), masks.begin(), masks.end());
    }

    /** Add a new bunch of masks - move semantics
     *
     * @tparam mask_index index for this mask type in masks container
     * @tparam mask_type mask type to be updated
     *
     * @param masks Vector of masks to be added
     *
     * @note in general can throw an exception
     */
    template<unsigned int mask_index,
             typename mask_type>
    void add_masks(vector_type<mask_type> &masks) noexcept(false) {
        // Get the mask group that will be updated
        auto &mask_group = std::get<mask_index>(_mask_tuple);
        // Reserve memory and copy new masks
        mask_group.reserve(mask_group.size() + masks.size());
        mask_group.insert(mask_group.end(), std::make_move_iterator(masks.begin()), std::make_move_iterator(masks.end()));
    }

    /** Append a mask store to the current one
     *
     * @param other The other mask store, move semantics
     *
     * @note in general can throw an exception
     */
    void append_masks(mask_store<vector_type, tuple_type, mask_types...> &&mask_store) noexcept(false) {
        return;
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
    inline void append_masks(mask_store<vector_type, tuple_type, mask_types...> &other) {
        // Add masks to current group
        auto &mask_group = std::get<current_index>(other);
        add_masks<current_index>(mask_group);

        // Next mask type
        if constexpr (current_index < std::tuple_size_v<mask_tuple> - 1) {
            return append_masks<current_index + 1>(other);
        }
    }

    private:

    /** tuple of masks */
    mask_tuple _mask_tuple;
};

}  // namespace detray
