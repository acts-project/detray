/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <tuple>
#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** @brief Struct that increments a given iterator in lockstep with the
 *        iteration index for convenience.
 *
 * @tparam container_t the container on which to perform the iteration
 * @tparam volume_t the type of volume from which to get the range in the
 *         container
 */
template <typename container_t,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_t>())),
          typename = decltype(std::end(std::declval<container_t>()))>
struct enumerator {
    /** Start position of iterator */
    size_t i;
    /** Iterator on container */
    container_type_iter iter;

    /** Build from start index and corresponding iterator - rvalue */
    DETRAY_HOST_DEVICE
    enumerator(size_t start, const container_type_iter &&iterator)
        : i(start), iter(iterator) {}

    /** Build from start index and corresponding iterator - lvalue */
    DETRAY_HOST_DEVICE
    enumerator(size_t start, const container_type_iter &iterator)
        : i(start), iter(iterator) {}

    /** Determine end of iteration */
    DETRAY_HOST_DEVICE
    bool operator!=(const enumerator &rhs) const { return iter != rhs.iter; }

    /** Increase index and iterator at once */
    DETRAY_HOST_DEVICE
    void operator++() {
        ++i;
        ++iter;
    }

    /** Tie them together for returning */
    DETRAY_HOST_DEVICE
    auto operator*() const { return std::tie(i, *iter); }
};

/** Helper utility to allow indexed enumeration with structured binding
 *
 * Usage:
 *
 * for (auto [ i, value ] = enumerate(container) ) { ... };
 *
 * with 'container' any stl-like container
 */
template <typename container_t,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_t>())),
          typename = decltype(std::end(std::declval<container_t>()))>
DETRAY_HOST_DEVICE constexpr auto enumerate(container_t &&iterable) {

    struct iterable_wrapper {
        container_t iterable;
        DETRAY_HOST_DEVICE
        decltype(auto) begin() {
            return enumerator<container_t>(0, std::begin(iterable));
        }
        DETRAY_HOST_DEVICE
        decltype(auto) end() {
            return enumerator<container_t>(0, std::end(iterable));
        }
    };

    return iterable_wrapper{std::forward<container_t>(iterable)};
}

}  // namespace detray
