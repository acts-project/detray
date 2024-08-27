/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

template <typename... Ts>
struct tuple {};

template <typename T, typename... Ts>
struct tuple<T, Ts...> {

    DETRAY_HOST_DEVICE constexpr tuple(){};

    template <
        typename U, typename... Us,
        std::enable_if_t<std::is_constructible_v<T, U &&> &&
                             std::is_constructible_v<tuple<Ts...>, Us &&...>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr tuple(const tuple<U, Us...> &o)
        : v(o.v), r(o.r) {}

    template <
        typename U, typename... Us,
        std::enable_if_t<std::is_constructible_v<T, U &&> &&
                             std::is_constructible_v<tuple<Ts...>, Us &&...>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr tuple(U &&_v, Us &&... _r)
        : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

    T v;
    tuple<Ts...> r;
};

template <std::size_t I, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const detray::tuple<Ts...> &t) noexcept {
    static_assert(I < sizeof...(Ts),
                  "Attempt to access index greater than tuple size.");

    if constexpr (I == 0) {
        return t.v;
    } else {
        return ::detray::get<I - 1>(t.r);
    }
}

template <std::size_t I, typename... Ts>
DETRAY_HOST_DEVICE auto &get(detray::tuple<Ts...> &t) noexcept {
    static_assert(I < sizeof...(Ts),
                  "Attempt to access index greater than tuple size.");

    if constexpr (I == 0) {
        return t.v;
    } else {
        return ::detray::get<I - 1>(t.r);
    }
}

template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const detray::tuple<T, Ts...> &t) noexcept {

    if constexpr (std::is_same_v<U, T>) {
        return t.v;
    } else if constexpr (sizeof...(Ts) > 0) {
        static_assert((std::is_same_v<U, Ts> || ...),
                      "Type not found in tuple.");
        return ::detray::get<U, Ts...>(t.r);
    }
}

template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE auto &get(detray::tuple<T, Ts...> &t) noexcept {

    if constexpr (std::is_same_v<U, T>) {
        return t.v;
    } else if constexpr (sizeof...(Ts) > 0) {
        static_assert((std::is_same_v<U, Ts> || ...),
                      "Type not found in tuple.");
        return ::detray::get<U, Ts...>(t.r);
    }
}

template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<Ts &...> tie(Ts &... args) {
    return detray::tuple<Ts &...>(args...);
}
}  // namespace detray
