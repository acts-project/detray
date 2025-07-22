/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#ifndef DETRAY_COMPILE_VITIS
#include <vecmem/memory/memory_resource.hpp>
#endif // DETRAY_COMPILE_VITIS

// System include(s)
#include <memory>
#include <type_traits>

namespace detray::detail {

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
    using buffer_type = dmulti_buffer<get_buffer_t<Ts>...>;

    /// Empty container - default alloc
    constexpr tuple_container() = default;
    /// Move constructor
    constexpr tuple_container(tuple_container &&) = default;
    /// Move assignment operator
    constexpr tuple_container &operator=(tuple_container &&) = default;

    /// Copy construct from element types
    constexpr explicit tuple_container(const Ts &... args) : _tuple(args...) {}

    /// Construct with a specific vecmem memory resource @param resource
    /// (host-side only)
#ifndef DETRAY_COMPILE_VITIS
    template <typename allocator_t = vecmem::memory_resource,
              std::enable_if_t<not is_device_view<allocator_t>::value , bool> = true>
    DETRAY_HOST explicit tuple_container(allocator_t &resource)
        : _tuple(Ts(&resource)...) {}
#endif // DETRAY_COMPILE_VITIS

    /// Copy Construct with a specific (vecmem) memory resource @param resource
    /// (host-side only)
#ifndef DETRAY_COMPILE_VITIS
    template <
        typename allocator_t = vecmem::memory_resource,
        typename T = tuple_t<Ts...>,
        std::enable_if_t<std::is_same<T, std::tuple<Ts...>>::value , bool> = true>
    DETRAY_HOST explicit tuple_container(allocator_t &resource,
                                         const Ts &... args)
        : _tuple(std::allocator_arg, resource, args...) {}
#endif // DETRAY_COMPILE_VITIS

    /// Construct from the container @param view type. Mainly used device-side.
    template <typename tuple_view_t,
              std::enable_if_t<is_device_view<tuple_view_t>::value , bool> = true>
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

    /// @returns the tuple element corresponding to the type @tparam T
    template <typename T>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() const noexcept {
        return detail::get<T>(_tuple);
    }

    /// @returns the tuple element corresponding to the type @tparam T
    template <typename T>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get() noexcept {
        return detail::get<T>(_tuple);
    }

    /// @returns a tuple of the views of all elements - non-const
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST auto get_data() -> view_type {
        return get_data(std::make_index_sequence<sizeof...(Ts)>{});
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns a tuple of the views of all elements - const
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST auto get_data() const -> const_view_type {
        return get_data(std::make_index_sequence<sizeof...(Ts)>{});
    }
#endif // DETRAY_COMPILE_VITIS

    /// Calls a functor with a all elements as parameters.
    ///
    /// @returns the functor result.
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE decltype(auto) apply(Args &&... As) const {

        return apply_impl<functor_t>(std::make_index_sequence<sizeof...(Ts)>{},
                                     std::forward<Args>(As)...);
    }

    /// Visits a tuple element according to its @param idx and calls
    /// @tparam functor_t with the arguments @param As on it.
    ///
    /// @returns the functor result (this is necessarily always of the same
    /// type, regardless the input tuple element type).
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE decltype(auto) visit(const std::size_t idx,
                                            Args &&... As) const {

        return visit<functor_t>(idx, std::make_index_sequence<sizeof...(Ts)>{},
                                std::forward<Args>(As)...);
    }

    private:
    /// @returns the view for all contained types.
#ifndef DETRAY_COMPILE_VITIS
    template <bool all_viewable = std::conjunction<detail::has_view<Ts>...>::value ,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST view_type get_data(std::index_sequence<I...> /*seq*/) noexcept {
        return view_type{detray::get_data(detail::get<I>(_tuple))...};
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns the const view for all contained types.
#ifndef DETRAY_COMPILE_VITIS
    template <bool all_viewable = std::conjunction<detail::has_view<Ts>...>::value ,
              std::size_t... I, std::enable_if_t<all_viewable, bool> = true>
    DETRAY_HOST const_view_type
    get_data(std::index_sequence<I...> /*seq*/) const noexcept {
        return const_view_type{detray::get_data(detail::get<I>(_tuple))...};
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns a tuple constructed from the elements @param view s.
    template <typename tuple_view_t, std::size_t... I,
              std::enable_if_t<is_device_view<tuple_view_t>::value , bool> = true>
    DETRAY_HOST_DEVICE auto unroll_views(tuple_view_t &view,
                                         std::index_sequence<I...> /*seq*/) {
        return detail::make_tuple<tuple_t>(Ts(detail::get<I>(view.m_view))...);
    }

    /// Variadic unrolling of the tuple that calls a functor on all elements of
    /// the tuple.
    ///
    /// @tparam functor_t functor that will be called on the elements.
    /// @tparam Args argument types for the functor
    template <typename functor_t, std::size_t... I, typename... Args>
    DETRAY_HOST_DEVICE decltype(auto) apply_impl(
        std::index_sequence<I...> /*seq*/, Args &&... As) const {

        // Call the functor on the tuple elements
        return functor_t{}(std::forward<Args>(As)..., get<I>()...);
    }

    /// Variadic unrolling of the tuple that calls a functor on the element that
    /// corresponds to @param idx.
    ///
    /// @tparam functor_t functor that will be called on the element.
    /// @tparam Args argument types for the functor
    /// @tparam first_idx Current index into the container tuple. Is converted
    ///         to an id_t and tested aginst the given id.
    /// @tparam remaining_idcs te remaining tuple indices to be tested.
    ///
    /// @see https://godbolt.org/z/qd6xns7KG
    template <typename functor_t, typename... Args, std::size_t first_idx,
              std::size_t... remaining_idcs>
    DETRAY_HOST_DEVICE std::invoke_result_t<
        functor_t, const detail::tuple_element_t<0, tuple_type> &, Args...>
    visit(const std::size_t idx,
          std::index_sequence<first_idx, remaining_idcs...> /*seq*/,
          Args &&... As) const {

        // Check if the first tuple index is matched to the target ID
        if (idx == first_idx) {
            return functor_t()(get<first_idx>(), std::forward<Args>(As)...);
        }
        // Check the next ID
        if constexpr (sizeof...(remaining_idcs) >= 1u) {
            return visit<functor_t>(idx,
                                    std::index_sequence<remaining_idcs...>{},
                                    std::forward<Args>(As)...);
        }
        // If there is no matching ID, return default output
        if constexpr (not std::is_same<
                          std::invoke_result_t<
                              functor_t,
                              const detail::tuple_element_t<0, tuple_type> &,
                              Args...>,
                          void>::value ) {
            return {};
        }
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

/// Get @c tuple_container element type
template <std::size_t N, template <typename...> class tuple_t, typename... Ts>
struct tuple_element<N, detail::tuple_container<tuple_t, Ts...>>
    : public detail::tuple_element<N, tuple_t<Ts...>> {};

/// Get @c tuple_container size
template <template <typename...> class tuple_t, typename... Ts>
struct tuple_size<detail::tuple_container<tuple_t, Ts...>>
    : public detail::tuple_size<tuple_t<Ts...>> {};

}  // namespace detray::detail
