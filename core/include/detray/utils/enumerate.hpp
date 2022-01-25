/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <tuple>
#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/indexing.hpp"
#include "detray/definitions/detail/accessor.hpp"

namespace detray {

/** @brief Struct that implements a range by providing start and end iterators.
 *
 * @tparam container_type the container on which to perform the iteration
 * @tparam volume_type the type of volume from which to get the range in the
 *         container
 */
template <typename container_type>
struct iterator_range {
    using container_type_iter =
        decltype(std::begin(std::declval<container_type>()));

    /** Delete default constructor */
    iterator_range() = delete;

    /** Always construct from a container and a range
     *
     * @param iterable container to iterate over
     * @param range start and end position for iteration
     */
    template <typename range_type>
    DETRAY_HOST_DEVICE iterator_range(const container_type &iterable,
                                      const range_type &range)
        : _start(std::begin(iterable) + detail::get<0>(range)),
          _end(std::begin(iterable) + detail::get<1>(range)) {}

    /** @return start position of range on container. */
    DETRAY_HOST_DEVICE
    inline auto &begin() { return _start; }

    /** @return end position of range on container. */
    DETRAY_HOST_DEVICE
    inline auto &end() { return _end; }

    /** Does this describe the same range? */
    DETRAY_HOST_DEVICE
    bool operator!=(const iterator_range &rhs) {
        return _start != rhs._start or _end != rhs._end;
    }

    /** @return element at position i, relative to iterator range. */
    DETRAY_HOST_DEVICE
    inline decltype(auto) operator[](const dindex i) { return *(_start + i); }

    /** @return element at position i, relative to iterator range - const */
    DETRAY_HOST_DEVICE
    inline decltype(auto) operator[](const dindex i) const {
        return *(_start + i);
    }

    /** @return the offset of the range start into the container. */
    DETRAY_HOST_DEVICE
    inline const auto offset(const container_type &iterable) {
        return std::distance(_start, std::begin(iterable));
    }

    /** Start and end position of a range */
    container_type_iter _start, _end;
};

/** @brief Struct that increments a given iterator in lockstep with the
 *        iteration index for convenience.
 *
 * @tparam container_type the container on which to perform the iteration
 * @tparam volume_type the type of volume from which to get the range in the
 *         container
 */
template <typename container_type,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_type>())),
          typename = decltype(std::end(std::declval<container_type>()))>
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

/** Get an interator range for the constituents of a volume on their container.
 *
 * @tparam container_type the container on which to perform the iteration
 * @tparam volume_type the type of volume from which to get the range in the
 *         container
 *
 * @param iterable reference to the container
 * @param volume the volume
 *
 * @returns an iterator range on the container, according to the volumes range
 */
template <
    typename container_type,
    typename = typename std::remove_reference_t<container_type>::value_type,
    typename volume_type,
    typename = typename std::remove_reference_t<volume_type>::volume_def>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_type &iterable,
                                               volume_type &&volume) {
    return iterator_range(
        iterable, volume.template range<typename container_type::value_type>());
}

/** Overload of the range-function for dindex_range */
template <typename container_type>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_type &iterable,
                                               const dindex_range &range) {

    return iterator_range(iterable, range);
}

/** Overload of the range-function for a single index */
template <typename container_type>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_type &iterable,
                                               const dindex &i) {

    return iterator_range(iterable, dindex_range{i, i + 1});
}

/** Helper utility to allow indexed enumeration with structured binding
 *
 * Usage:
 *
 * for (auto [ i, value ] = enumerate(container) ) { ... };
 *
 * with 'container' any stl-like container
 */
template <typename container_type,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_type>())),
          typename = decltype(std::end(std::declval<container_type>()))>
DETRAY_HOST_DEVICE constexpr auto enumerate(container_type &&iterable) {

    struct iterable_wrapper {
        container_type iterable;
        DETRAY_HOST_DEVICE
        decltype(auto) begin() {
            return enumerator<container_type>(0, std::begin(iterable));
        }
        DETRAY_HOST_DEVICE
        decltype(auto) end() {
            return enumerator<container_type>(0, std::end(iterable));
        }
    };

    return iterable_wrapper{std::forward<container_type>(iterable)};
}

/** Helper utility to allow indexed enumeration with structured binding for
 *  a given range
 *
 * Usage:
 *
 * for (auto [ i, value ] = enumerate(container, range) ) { ... };
 *
 * with 'container' any stl-like container and r a range type that is compatible
 * with std::get or a type for which an overload to the range() function exists.
 */
template <typename container_type, typename range_type,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_type>())),
          typename = decltype(std::end(std::declval<container_type>()))>
DETRAY_HOST_DEVICE constexpr inline auto enumerate(
    const container_type &iterable, range_type &&r) {

    struct iterable_wrapper {
        const container_type_iter iter;
        iterator_range<container_type> range_iter;

        DETRAY_HOST_DEVICE
        decltype(auto) begin() {
            return enumerator<container_type>(
                static_cast<size_t>(std::distance(iter, range_iter.begin())),
                range_iter.begin());
        }
        DETRAY_HOST_DEVICE
        decltype(auto) end() {
            return enumerator<container_type>(
                static_cast<size_t>(std::distance(iter, range_iter.begin())),
                range_iter.end());
        }
    };

    return iterable_wrapper{std::begin(iterable),
                            range(iterable, std::forward<range_type>(r))};
}

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
    typename range_type,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<0>(std::declval<range_type>()))>>>,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<1>(std::declval<range_type>()))>>>>
DETRAY_HOST_DEVICE constexpr auto sequence(range_type range) {

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
        range_type _iterable;

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

    return iterable_wrapper{std::forward<range_type>(range)};
}

}  // namespace detray
