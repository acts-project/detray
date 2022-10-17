/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <type_traits>

namespace detray {

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
    template <typename allocator_t = vecmem::memory_resource,
              std::enable_if_t<not is_device_view_v<allocator_t>, bool> = true>
    DETRAY_HOST explicit tuple_container(allocator_t &resource)
        : _tuple(Ts(&resource)...) {}

    /// Copy Construct with a specific (vecmem) memory resource @param resource
    /// (host-side only)
    template <
        typename allocator_t = vecmem::memory_resource,
        typename T = tuple_t<Ts...>,
        std::enable_if_t<std::is_same_v<T, std::tuple<Ts...>>, bool> = true>
    DETRAY_HOST explicit tuple_container(allocator_t &resource,
                                         const Ts &...args)
        : _tuple(std::allocator_arg, resource, args...) {}

    /// Construct from the container @param view type. Mainly used device-side.
    template <typename tuple_view_t,
              std::enable_if_t<is_device_view_v<tuple_view_t>, bool> = true>
    DETRAY_HOST_DEVICE tuple_container(tuple_view_t &view)
        : _tuple(
              unroll_views(view, std::make_index_sequence<sizeof...(Ts)>{})) {}

    /// @returns the size of the tuple
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        return sizeof...(Ts);
    }

    /// @returns the tuple element corresponding to the index @tparam idx
    template <std::size_t idx>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() const noexcept {
        return detail::get<idx>(_tuple);
    }

    /// @returns the tuple element corresponding to the index @tparam idx
    template <std::size_t idx>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() noexcept {
        return detail::get<idx>(_tuple);
    }

    /// @returns the tuple element corresponding to the index @tparam idx
    template <typename T>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() const noexcept {
        return detail::get<T>(_tuple);
    }

    /// @returns the tuple element corresponding to the index @tparam idx
    template <typename T>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() noexcept {
        return detail::get<T>(_tuple);
    }

    /// @returns a tuple of the views of all elements - non-const
    DETRAY_HOST auto get_data() -> view_type {
        return get_data(std::make_index_sequence<sizeof...(Ts)>{});
    }

    /// @returns a tuple of the views of all elements - const
    DETRAY_HOST auto get_data() const -> const_view_type {
        return get_data(std::make_index_sequence<sizeof...(Ts)>{});
    }

    /// Calls a functor with an element with a specific index.
    ///
    /// @return the functor output
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type call(
        const std::size_t idx, Args &&...As) const {

        return unroll_call<functor_t>(idx,
                                      std::make_index_sequence<sizeof...(Ts)>{},
                                      std::forward<Args>(As)...);
    }

    private:
    /// @returns the view for all contained types.
    template <bool all_viewable = std::conjunction_v<detail::get_view<Ts>...>,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST view_type get_data(std::index_sequence<I...> /*seq*/) noexcept {
        return {detray::get_data(detail::get<I>(_tuple))...};
    }

    /// @returns the const view for all contained types.
    template <bool all_viewable = std::conjunction_v<detail::get_view<Ts>...>,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST const_view_type
    get_data(std::index_sequence<I...> /*seq*/) const noexcept {
        return {detray::get_data(detail::get<I>(_tuple))...};
    }

    /// @returns a tuple constructed from the elements @param view s.
    template <typename tuple_view_t, std::size_t... I,
              std::enable_if_t<is_device_view_v<tuple_view_t>, bool> = true>
    DETRAY_HOST_DEVICE auto unroll_views(tuple_view_t &view,
                                         std::index_sequence<I...> /*seq*/) {
        return detail::make_tuple<tuple_t>(Ts(detail::get<I>(view.m_view))...);
    }

    /// Variadic unrolling of the tuple that calls a functor with the element.
    ///
    /// @tparam functor_t functor that will be called on the element.
    /// @tparam Args argument types for the functor
    /// @tparam first_idx Current index into the container tuple. Is converted
    ///         to an id_t and tested aginst the given id.
    /// @tparam remaining_idcs te remaining tuple indices to be tested.
    template <typename functor_t, typename... Args, std::size_t first_idx,
              std::size_t... remaining_idcs>
    DETRAY_HOST_DEVICE typename functor_t::output_type unroll_call(
        const std::size_t idx,
        std::index_sequence<first_idx, remaining_idcs...> /*seq*/,
        Args &&...As) const {

        // Check if the first tuple index is matched to the target ID
        if (idx == first_idx) {
            const auto &elem = get<first_idx>();

            return functor_t()(elem, std::forward<Args>(As)...);
        }
        // Check the next ID
        if constexpr (sizeof...(remaining_idcs) >= 1) {
            return unroll_call<functor_t>(
                idx, std::index_sequence<remaining_idcs...>{},
                std::forward<Args>(As)...);
        }
        // If there is no matching ID, return null output
        return typename functor_t::output_type{};
    }

    /// The underlying tuple container
    tuple_type _tuple;
};

/// Overloads to 'get' for the tuple container
/// @{
template <std::size_t idx, template <typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const detail::tuple_container<tuple_t, Ts...> &container) {
    return container.template get<idx>();
}

template <std::size_t idx, template <typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    detail::tuple_container<tuple_t, Ts...> &container) {
    return container.template get<idx>();
}

template <typename T, template <typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const detail::tuple_container<tuple_t, Ts...> &container) {
    return container.template get<T>();
}

template <typename T, template <typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    detail::tuple_container<tuple_t, Ts...> &container) {
    return container.template get<T>();
}
/// @}

}  // namespace detail

}  // namespace detray