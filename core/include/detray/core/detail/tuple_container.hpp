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
        return detail::get<ID>(m_container).size();
    }

    /**
     * Return if the tuple element is empty.
     *
     * @tparam ID is the index of tuple element
     * @return true if the tuple element is empty
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE bool empty() const {
        return detail::get<ID>(m_container).empty();
    }

    /**
     * Return a tuple element (non-const access)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE constexpr auto &group() {
        return detail::get<ID>(m_container);
    }

    /**
     * Return a tuple element (const access)
     *
     * @tparam mask_id index of requested mask type in masks container
     * @return vector of masks of a given type.
     */
    template <id_t ID>
    DETRAY_HOST_DEVICE constexpr const auto &group() const {
        return detail::get<ID>(m_container);
    }

    /** Enforce usage id_t in the code and do some (limited)
     *  checking.
     *
     * @tparam ref_idx matches to index arg to perform static checks
     * @param index argument to be converted to valid id type
     *
     * @return the matching ID type.
     */
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr id_t to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                ref_idx < sizeof...(Ts),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<id_t>(index);
        }
        if constexpr (ref_idx < sizeof...(Ts) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<id_t>(sizeof...(Ts));
    }

    /** Execute functor for a group with specific ID. The group is found by
     * unrolling varidically
     *
     * @tparam functor_t is the functor type
     * @tparam size_type is type for index
     * @tparam Args is argument type for the functor
     *
     * @param id is the target group index
     * @param As is the functor arguments
     *
     * @return the functor output
     */
    template <typename functor_t, typename size_type, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type execute(
        const size_type id, Args &&... As) const {

        return unroll<functor_t>(id, std::make_index_sequence<sizeof...(Ts)>{},
                                 std::forward<Args>(As)...);
    }

    protected:
    container_type m_container;

    private:
    /** Variadic unroll function used for execute function
     */
    template <typename functor_t, typename size_type, typename... Args,
              std::size_t first_id, std::size_t... remaining_ids>
    DETRAY_HOST_DEVICE typename functor_t::output_type unroll(
        const size_type id,
        std::index_sequence<first_id, remaining_ids...> /*seq*/,
        Args &&... As) const {

        // Check if the first ID is matched to the target ID
        if (id == first_id) {
            const auto &gr = this->group<to_id(first_id)>();

            return functor_t()(gr, std::forward<Args>(As)...);
        }

        // Check the next ID
        if constexpr (sizeof...(remaining_ids) >= 1) {
            return unroll<functor_t>(id,
                                     std::index_sequence<remaining_ids...>{},
                                     std::forward<Args>(As)...);
        }

        // If there is no matching ID, return null output
        return typename functor_t::output_type{};
    }
};

}  // namespace detray
