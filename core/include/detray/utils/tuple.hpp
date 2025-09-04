/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>
#include <utility>

#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {
template <typename... Ts>
struct tuple {
    DETRAY_HOST_DEVICE constexpr tuple() {}
};

template <typename T, typename... Ts>
struct tuple<T, Ts...> {
    DETRAY_HOST_DEVICE constexpr tuple() {}

    template <typename U, typename... Us>
    DETRAY_HOST_DEVICE constexpr tuple(const tuple<U, Us...> &o)
        : v(o.v), r(o.r) {}

    template <
        typename U, typename... Us,
        std::enable_if_t<std::is_constructible<T, U &&>::value  &&
                             std::is_constructible<tuple<Ts...>, Us &&...>::value ,
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

template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<Ts &...> tie(Ts &... args) {
    return detray::tuple<Ts &...>(args...);
}
}  // namespace detray
