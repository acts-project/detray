/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_traits.hpp"
#include "detray/core/detail/vecmem_types.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

// System include(s)
#include <memory>
#include <type_traits>

namespace detray::detail {

template<typename Test, template<typename...> class Ref>
struct is_specialization : std::false_type {};

template<template<typename...> class Ref, typename... Args>
struct is_specialization<Ref<Args...>, Ref>: std::true_type {};

/// @tparam Ts types of tuple elements.
template <typename... Ts>
struct tuple_host {
    DETRAY_HOST constexpr tuple_host() {}
};

template <typename T, typename... Ts>
struct tuple_host<T, Ts...> {
    using h_t = typename T::host;

    DETRAY_HOST constexpr tuple_host() {}

    template <typename U, typename... Us>
    DETRAY_HOST constexpr tuple_host(const tuple<U, Us...> &o)
        : v(o.v), r(o.r) {}

    template <
        typename U, typename... Us
        ,std::enable_if_t<std::is_constructible_v<h_t, U &&> &&
                             std::is_constructible_v<tuple_host<Ts...>, Us &&...>,
                         bool> = true
                         >
    DETRAY_HOST constexpr tuple_host(U &&_v, Us &&... _r)
        : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

    h_t v;
    tuple_host<Ts...> r;
};

/// @tparam Ts types of tuple elements.
template <typename T, typename... Ts>
struct tuple_buffer {
    using b_t = typename T::buffer;

    template <typename size, typename... sizes>
    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource, size sz, sizes... szs)
        : v([&]{
            if constexpr(is_specialization<b_t, vecmem::data::vector_buffer>::value) {
                return b_t(sz, resource);
            }
            else {
                return b_t(vecmem::make_unique_alloc<typename b_t::element_type>(resource));
            }
        }()),
        r([&]{
            if constexpr(is_specialization<b_t, vecmem::data::vector_buffer>::value) {
                return tuple_buffer<Ts...>(resource, szs...);
            }
            else {
                return tuple_buffer<Ts...>(resource, sz, szs...);
            }
        }()) {}

    b_t v;
    tuple_buffer<Ts...> r;
};

/// @tparam T type of last tuple element.
template <typename T>
struct tuple_buffer<T> {
    using b_t = typename T::buffer;

    template <typename size>
    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource, size sz)
        : v(sz, resource) {
            static_assert(is_specialization<b_t, vecmem::data::vector_buffer>::value, "Number of given sizes larger than number of vectors.");
        }

    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource) 
        : v(vecmem::make_unique_alloc<typename b_t::element_type>(resource)) {
            static_assert(!is_specialization<b_t, vecmem::data::vector_buffer>::value, "Number of given sizes smaller than number of vectors.");
        }
    b_t v;
};

template <typename T, typename... Ts>
struct tuple_view {
    using v_t = typename T::view;

    DETRAY_HOST constexpr tuple_view(tuple_buffer<T, Ts...> &buff) :
        v(test::get_data(buff.v)), r(buff.r) {}

    v_t v;
    tuple_view<Ts...> r;
};

template <typename T>
struct tuple_view<T> {
    using v_t = typename T::view;

    DETRAY_HOST constexpr tuple_view(tuple_buffer<T> &buff) :
        v(test::get_data(buff.v)) {}

    v_t v;
};

template <typename T, typename... Ts>
struct tuple_device {
    using d_t = typename T::device;

    DETRAY_HOST constexpr tuple_device(tuple_view<T, Ts...> &view) :
        v(view.v), r(view.r) {}

    d_t v;
    tuple_device<Ts...> r;
};

template <typename T>
struct tuple_device<T> {
    using d_t = typename T::device;

    DETRAY_HOST constexpr tuple_device(tuple_view<T> &view) :
        v(view.v) {}

    d_t v;
};

/// Tuple types in a vecmem-like interface
/// Their capable of handling vecmem::vector like objects as well as POD.
///
/// Host type enables typical memory access, dynamic allocation, etc.
/// Buffer type allows memory allocation on different devices
/// View type allows passing buffers via non-ownership to CPU kernels
/// Device type allows reading/writing onto view types
template <typename... Ts>
struct tuple_types {
    using host = tuple_host<Ts...>;
    using buffer = tuple_buffer<Ts...>;
    using view = tuple_view<Ts...>;
    using device = tuple_device<Ts...>;
};

template <std::size_t I, template<typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const tuple_t<Ts...> &t) noexcept {
    static_assert(I < sizeof...(Ts),
                  "Attempt to access index greater than tuple size.");

    if constexpr (I == 0) {
        return t.v;
    } else {
        return get<I - 1>(t.r);
    }
}

template <std::size_t I, template<typename...> class tuple_t, typename... Ts>
DETRAY_HOST_DEVICE auto &get(tuple_t<Ts...> &t) noexcept {
    static_assert(I < sizeof...(Ts),
                  "Attempt to access index greater than tuple size.");

    if constexpr (I == 0) {
        return t.v;
    } else {
        return get<I - 1>(t.r);
    }
}

}  // namespace detray::detail
