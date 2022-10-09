/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray includes
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

namespace detail {

struct viewable_container_base {};

/// Interface definition for detray data containers that need to be stored in
/// containers and have to do their own host-device memory management.
/// @{
template <typename T, class, typename = void>
struct viewable_container;

/// The @tparam T does not perform memory management of its own, but should be
/// part of a container that does.
template <typename T, template <typename> class container_t>
    struct viewable_container < T,
    container_t, std::enable_if_t<not is_viewable_v<T>>,
    void >> : public viewable_container_base, public container_t<T> {

    /// Container types derived from the container given
    /// @{
    using value_type = T;
    using reference = T &;
    using iterator = typename container_t<T>::iterator;
    using const_iterator = typename container_t<T>::const_iterator;
    using difference_type =
        typename std::iterator_traits<iterator>::difference_type;
    using size_type = typename container_t<T>::size_type;
    /// @}

    using view_type = get_view_t<T>;
}

/// The @tparam container_t is already a container with a view type defined.
template <typename container_t>
struct viewable_container<
    container_t, void,
    std::enable_if_t<
        not std::is_base_of<viewable_container_base, container_t> and
            is_viewable_v<container_t>,
        void>> : public viewable_container_base,
                 public container_t {

    /// Container types derived from the container of objects
    /// @{
    using value_type = typename container_t::value_type;
    using reference = typename container_t::value_type &;
    using iterator = decltype(std::begin(std::declval<container_t &>()));
    using const_iterator = decltype(std::cbegin(std::declval<container_t &>()));
    using difference_type =
        typename std::iterator_traits<iterator>::difference_type;
    using size_type = typename container_t::size_type;
    /// @}

    using view_type = typename container_t::view_type;
}

/// @}

}  // namespace detail

/// Get the vecmem view object of a detray viewable data container
template <typename T>
get_data(viewable_container<T> coll &) {
    return get_data(coll);
}

/// Get the vecmem view object of a detray viewable data container
template <typename T>
get_data(const viewable_container<T> coll &) {
    return get_data(coll);
}

}  // namespace detray