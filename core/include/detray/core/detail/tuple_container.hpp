/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Tuple based container.
 *
 * @tparam tuple_t is the type of tuple
 * @tparam id_t is an enum that is compared to an indexing integer
 * @tparam Ts... is the types of tuple elements
 */
template <template <typename...> class tuple_t, typename id_t, typename... Ts>
class tuple_container {

    public:
    // Convenient type declarations
    static constexpr std::size_t m_tuple_size = sizeof...(Ts);
    template <typename... Args>
    using tuple_type = tuple_t<Args...>;
    using id_type = id_t;

    /**
     * Container type definitions.
     *
     * container_type does not follow the tuple_type.
     * vtuple, defined in accessor.hpp, has different types based on the file
     * location 1) std::tuple in *.cpp/hpp; 2) thrust::tuple in *.cu
     *
     * It is definitely desirable to use tuple_type but it causes illegal memory
     * access error in CUDA applications (See also Detray Issue #154.)
     */
    using container_type = vtuple::tuple<Ts...>;

    /**
     * Move constructor
     *
     * @param container is lvalue input container
     */
    DETRAY_HOST_DEVICE
    tuple_container(container_type &&container)
        : m_container(std::move(container)) {}

    /**
     * Get container (const access).
     *
     * @return the container
     */
    DETRAY_HOST_DEVICE
    const auto &get() const { return m_container; }

    /**
     * Get container (non-const access).
     *
     * @return the container
     */
    DETRAY_HOST_DEVICE
    auto &get() { return m_container; }

    /**
     * Get the size of tuple.
     *
     * @return the size of tuple
     */
    DETRAY_HOST_DEVICE constexpr std::size_t size() const {
        return m_tuple_size;
    }

    /**
     * Get the size of tuple element.
     *
     * @return the size of tuple element
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE size_t size() const {
        return detail::get<to_index(ID)>(m_container).size();
    }

    /**
     * Return if the tuple element is empty.
     *
     * @tparam ID is the index of tuple element
     * @return true if the tuple element is empty
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE bool empty() const {
        return detail::get<to_index(ID)>(m_container).empty();
    }

    /**
     * Return a tuple element (non-const access)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE constexpr auto &group() {
        return detail::get<to_index(ID)>(m_container);
    }

    /**
     * Return a tuple element (const access)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE constexpr const auto &group() const {
        return detail::get<to_index(ID)>(m_container);
    }

    /// @brief Convert index of the tuple to an id.
    ///
    /// The function uses unrolling, so that it can be used as a constant expr.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr id_t to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                ref_idx < sizeof...(Ts),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in the tuple container.");
            return static_cast<id_t>(index);
        }
        if constexpr (ref_idx < sizeof...(Ts) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<id_t>(sizeof...(Ts));
    }

    /// @brief Convert an id to an index of the tuple.
    ///
    /// The function uses unrolling, so that it can be used as a constant expr.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param id type id that should be used to index a tuple element
    ///
    /// @return the matching index.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr std::size_t to_index(const id_t id) {
        if (to_id(ref_idx) == id) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                ref_idx < sizeof...(Ts),
                "Index out of range: This ID cannot be used to index the tuple "
                "container.");
            return ref_idx;
        }
        if constexpr (ref_idx < sizeof...(Ts) - 1) {
            return to_index<ref_idx + 1>(id);
        }
        // This produces a compiler error when used in type unrolling code
        return sizeof...(Ts);
    }

    /// Calls a functor for a group with specific ID. The group is found by
    /// unrolling varidically
    ///
    /// @tparam functor_t functor that will be called on the group.
    /// @tparam Args argument types for the functor
    ///
    /// @param id is the target group id
    /// @param As additional functor arguments
    ///
    /// @return the functor output
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type call(
        const id_t id, Args &&...As) const {

        // An invalid range will be interpreted by the detray range iterator to
        // mean the entire range. Otherwise use overload function below to
        // specify a valid range
        return unroll<functor_t>(id, dindex_range{0, dindex_invalid},
                                 std::make_index_sequence<sizeof...(Ts)>{},
                                 std::forward<Args>(As)...);
    }

    /// Calls a functor for a group with specific ID. The group is found by
    /// unrolling varidically
    ///
    /// @tparam functor_t functor that will be called on the group.
    /// @tparam link_t how to reference a group and its entries.
    /// @tparam Args argument types for the functor
    ///
    /// @param link contains the group id and an index into the group
    /// @param As additional functor arguments
    ///
    /// @return the functor output
    template <typename functor_t, typename link_t, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type call(
        const link_t link, Args &&...As) const {

        return unroll<functor_t>(detail::get<0>(link), detail::get<1>(link),
                                 std::make_index_sequence<sizeof...(Ts)>{},
                                 std::forward<Args>(As)...);
    }

    protected:
    container_type m_container;

    private:
    /// Variadic unroll function used for execute function
    ///
    /// @tparam functor_t functor that will be called on the group (members).
    /// @tparam index_t how to reference a member(s) in the group. Can be a
    ///         single index/range/multiindex
    /// @tparam Args argument types for the functor
    /// @tparam first_idx Current index into the container tuple. Is converted
    ///         to an id_t and tested aginst the given id.
    /// @tparam remaining_idcs te remaining tuple indices to be tested.
    template <typename functor_t, typename index_t, typename... Args,
              std::size_t first_idx, std::size_t... remaining_idcs>
    DETRAY_HOST_DEVICE typename functor_t::output_type unroll(
        const id_t id, const index_t index,
        std::index_sequence<first_idx, remaining_idcs...> /*seq*/,
        Args &&...As) const {

        // Check if the first tuple index is matched to the target ID
        if (id == to_id(first_idx)) {
            const auto &gr = this->group<to_id(first_idx)>();

            return functor_t()(gr, index, std::forward<Args>(As)...);
        }

        // Check the next ID
        if constexpr (sizeof...(remaining_idcs) >= 1) {
            return unroll<functor_t>(id, index,
                                     std::index_sequence<remaining_idcs...>{},
                                     std::forward<Args>(As)...);
        }

        // If there is no matching ID, return null output
        return typename functor_t::output_type{};
    }
};

}  // namespace detray
