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
    DETRAY_HOST constexpr tuple_host() {}

    template <typename U, typename... Us>
    DETRAY_HOST constexpr tuple_host(const tuple<U, Us...> &o)
        : v(o.v), r(o.r) {}

    template <
        typename U, typename... Us,
        std::enable_if_t<std::is_constructible_v<T, U &&> &&
                             std::is_constructible_v<tuple_host<Ts...>, Us &&...>,
                         bool> = true>
    DETRAY_HOST constexpr tuple_host(U &&_v, Us &&... _r)
        : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

    T v;
    tuple_host<Ts...> r;
};

/// @tparam Ts types of tuple elements.
template <typename T, typename... Ts>
struct tuple_buffer {
    template <typename size, typename... sizes>
    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource, size sz, sizes... szs)
        : v([&]{
            // This fails for the unique_alloc_ptr
            // static_assert(is_specialization<T, vecmem::data::vector_buffer>::value 
            //     || is_specialization<T, vecmem::unique_alloc_ptr>::value);
            if constexpr(is_specialization<T, vecmem::data::vector_buffer>::value) {
                return T(sz, resource);
            }
            else {
                return T(vecmem::make_unique_alloc<typename T::element_type>(resource));
            }
        }()),
        r([&]{
            if constexpr(is_specialization<T, vecmem::data::vector_buffer>::value) {
                return tuple_buffer<Ts...>(resource, szs...);
            }
            else {
                return tuple_buffer<Ts...>(resource, sz, szs...);
            }
        }()) {}

    T v;
    tuple_buffer<Ts...> r;
};

/// @tparam T type of last tuple element.
template <typename T>
struct tuple_buffer<T> {
    template <typename size>
    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource, size sz)
        : v(sz, resource) {
            static_assert(is_specialization<T, vecmem::data::vector_buffer>::value, "Number of given sizes larger than number of vectors.");
        }

    DETRAY_HOST constexpr tuple_buffer(vecmem::memory_resource& resource) 
        : v(resource) {
            static_assert(!is_specialization<T, vecmem::data::vector_buffer>::value, "Number of given sizes smaller than number of vectors.");
        }
    T v;
};

}  // namespace detray::detail
