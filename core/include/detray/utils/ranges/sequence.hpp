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

/** Helper method to (fake a) run over a single entry
 *
 * Usage:
 * for (auto i : sequence(j)) {}
 *
 * with j an unsinged dindex type
 *
 * @note sequence(2) will produce {2}
 *
 **/
DETRAY_HOST_DEVICE
constexpr auto sequence(dindex iterable) {

    struct iterator {
        /** Start and end of sequence */
        dindex i;
        size_t end;

        /** Determine whether we reach end of sequence */
        DETRAY_HOST_DEVICE
        bool operator!=(const iterator &rhs) const { return i != rhs.end; }

        /** Increase index and iterator at once */
        DETRAY_HOST_DEVICE
        void operator++() { ++i; }

        /** Tie them together for returning */
        DETRAY_HOST_DEVICE
        auto operator*() const { return i; }
    };

    /** Wrap up for iteration */
    struct iterable_wrapper {
        dindex iterable;

        DETRAY_HOST_DEVICE
        auto begin() { return iterator{iterable, iterable + 1}; }
        DETRAY_HOST_DEVICE
        auto end() { return iterator{iterable + 1, iterable + 1}; }
    };

    return iterable_wrapper{std::forward<dindex>(iterable)};
}

/** Helper method to run over a range
 *
 * Usage:
 * for (auto i : sequence(r)) {}
 *
 * with r an index range (integral) that can be accessed with std::get
 *
 * @note sequence({2,4}) will produce { 2, 3, 4 }
 *
 **/
template <
    typename range_t,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<0>(std::declval<range_t>()))>>>,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<1>(std::declval<range_t>()))>>>>
DETRAY_HOST_DEVICE constexpr auto sequence(range_t range) {

    struct iterator {
        /** Start and end of sequence */
        size_t i;
        size_t end;

        /** Determine whether we reach end of sequence */
        DETRAY_HOST_DEVICE
        bool operator!=(const iterator &rhs) const { return i != rhs.end; }

        /** Increase index and iterator at once */
        DETRAY_HOST_DEVICE
        void operator++() { ++i; }

        /** Tie them together for returning */
        DETRAY_HOST_DEVICE
        auto operator*() const { return i; }
    };

    /** Wrap up for iteration */
    struct iterable_wrapper {
        range_t _iterable;

        DETRAY_HOST_DEVICE
        auto begin() {
            return iterator{std::get<0>(_iterable), std::get<1>(_iterable)};
        }
        DETRAY_HOST_DEVICE
        auto end() {
            return iterator{std::get<1>(_iterable) + 1,
                            std::get<1>(_iterable) + 1};
        }
    };

    return iterable_wrapper{std::forward<range_t>(range)};
}

}  // namespace detray
