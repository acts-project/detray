/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <type_traits>

namespace detray {

using vecmem::get_data;

namespace detail {

/// @brief detray tuple wrapper.
///
/// @tparam An enum of type IDs that needs to match the [value] types of the
/// @c Ts pack.
/// @tparam tuple_t is the type of the underlying tuple container
/// @tparam Ts are the types of tuple elements. They need to define their own
/// vecmem view type in order to be managed by the container when moving data
/// from host to device
template <template <typename...> class tuple_t, typename... Ts>
class tuple_container {

    public:
    using tuple_type = tuple_t<detray::detail::unwrap_decay_t<Ts>...>;
    using view_type = dmulti_view<get_view_t<Ts>...>;
    using const_view_type = dmulti_view<get_view_t<const Ts>...>;

    /// Empty container - default alloc
    constexpr tuple_container() = default;

    /// Copy construct from element types
    constexpr explicit tuple_container(const Ts &...args) : _tuple(args...) {}

    /// Construct with a specific vecmem memory resource @param resource
    /// (host-side only)
    DETRAY_HOST explicit tuple_container(vecmem::memory_resource &resource)
        : _tuple(Ts(&resource)...) {}

    /// Copy Construct with a specific vecmem memory resource @param resource
    /// (host-side only)
    template <
        typename allocator_t = vecmem::memory_resource,
        typename T = tuple_t<Ts...>,
        std::enable_if_t<std::is_same_v<T, std::tuple<Ts...>>, bool> = true>
    DETRAY_HOST explicit tuple_container(allocator_t &resource,
                                         const Ts &...args)
        : _tuple(std::allocator_arg, resource, args...) {}

    /// Construct from the container view type. Mainly used device-side.
    ///
    /// @tparam is the type of input data container
    ///
    /// @param container_data is the data container
    template <typename tuple_view_t,
              std::enable_if_t<is_device_view_v<tuple_view_t>, bool> = true>
    DETRAY_HOST_DEVICE tuple_container(tuple_view_t &view)
        : _tuple(
              unroll_views(view, std::make_index_sequence<sizeof...(Ts)>{})) {}

    /// @returns the view for all contained types.
    template <bool all_viewable = std::conjunction_v<detail::get_view<Ts>...>,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST view_type get_data(std::index_sequence<I...> /*seq*/) {
        return {detray::get_data(detail::get<I>(_tuple))...};
    }

    /// @returns the view for all contained types.
    template <bool all_viewable = std::conjunction_v<detail::get_view<Ts>...>,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST const_view_type
    get_data(std::index_sequence<I...> /*seq*/) const {
        return {detray::get_data(detail::get<I>(_tuple))...};
    }

    /// @returns the size of a data collection by id
    DETRAY_HOST_DEVICE
    constexpr auto size() const -> std::size_t {
        return detail::tuple_size<tuple_type>::value;
    }

    private:
    /// Construct from the container view type. Mainly used device-side.
    ///
    /// @tparam is the type of input data container
    ///
    /// @param container_data is the data container
    template <typename tuple_view_t, std::size_t... I,
              std::enable_if_t<is_device_view_v<tuple_view_t>, bool> = true>
    DETRAY_HOST_DEVICE auto unroll_views(tuple_view_t &view,
                                         std::index_sequence<I...> /*seq*/) {
        return detail::make_tuple<tuple_t>(
            Ts(detail::get<I>(view.m_views).m_view)...);
    }

    tuple_type _tuple;
};

}  // namespace detail

/// A stand-alone function to get the vecmem view of the tuple container
///
/// @note the @c view_type typedef will not be available, if one of the element
/// types does not define a vecmem view.
///
/// @return the view on this tuple container
template <template <typename...> class tuple_t, typename... Ts>
inline typename detail::tuple_container<tuple_t, Ts...>::view_type get_data(
    detail::tuple_container<tuple_t, Ts...> &container) {
    return container.get_data(std::make_index_sequence<sizeof...(Ts)>{});
}

}  // namespace detray