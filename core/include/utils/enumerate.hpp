/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <tuple>
#include <type_traits>

#include "utils/indexing.hpp"

namespace detray {

/** @brief Struct that implements a range by providing start and end iterators.
 */
template <typename container_type,
          typename container_type_iter =
              decltype(std::cbegin(std::declval<container_type>())),
          typename = decltype(std::cend(std::declval<container_type>()))>
struct iterator_range {
    /** Delete default constructor */
    iterator_range() = delete;

    /** Always construct from a container and a range
     *
     * @param iterable container to iterate over
     * @param range start and end position for iteration
     */
    template <typename range_type>
    iterator_range(const container_type &iterable, range_type &&range)
        : _start(std::next(std::cbegin(iterable),
                           std::get<0>(std::forward<range_type>(range)))),
          _end(std::next(std::cbegin(iterable),
                         std::get<1>(std::forward<range_type>(range)))) {}

    /** @return start position of range on container. */
    inline auto &begin() { return _start; }

    /** @return end position of range on container. */
    inline auto &end() { return _end; }

    /** Does this describe the same range? */
    bool operator!=(const iterator_range &rhs) {
        return _start != rhs._start or _end != rhs._end;
    }

    /** @return element at position i, relative to iterator range. */
    inline decltype(auto) operator[](const dindex i) { return *(_start + i); }

    /** @return element at position i, relative to iterator range - const */
    inline decltype(auto) operator[](const dindex i) const {
        return *(_start + i);
    }

    inline const auto offset(const container_type &iterable) {
        return std::distance(_start, std::cbegin(iterable));
    }
    /** Start and end position of a range */
    container_type_iter _start, _end;
};

template <typename container_type, typename volume_type>
inline constexpr decltype(auto) range(const container_type &iterable,
                                      volume_type &&volume) {
    return iterator_range(
        iterable, volume.template range<typename container_type::value_type>());
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
constexpr auto enumerate(container_type &&iterable) {
    struct iterator {
        size_t i;
        container_type_iter iter;

        bool operator!=(const iterator &rhs) const { return iter != rhs.iter; }

        /** Increase index and iterator at once */
        void operator++() {
            ++i;
            ++iter;
        }

        /** Tie them together for returning */
        auto operator*() const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper {
        container_type iterable;
        auto begin() { return iterator{0, std::begin(iterable)}; }
        auto end() { return iterator{0, std::end(iterable)}; }
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
              decltype(std::cbegin(std::declval<container_type>())),
          typename = decltype(std::cend(std::declval<container_type>()))>
constexpr inline auto enumerate(const container_type &iterable,
                                range_type &&r) {

    struct iterator {
        size_t i;
        container_type_iter &iter;

        bool operator!=(const iterator &rhs) const { return iter != rhs.iter; }

        /** Increase index and iterator at once */
        void operator++() {
            ++i;
            ++iter;
        }

        /** Tie them together for returning */
        auto operator*() const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper {
        const container_type_iter iter;
        iterator_range<container_type> range_iter;

        auto begin() {
            return iterator{
                static_cast<size_t>(std::distance(iter, range_iter.begin())),
                range_iter.begin()};
        }
        auto end() {
            return iterator{
                static_cast<size_t>(std::distance(iter, range_iter.begin())),
                range_iter.end()};
        }
    };
    return iterable_wrapper{std::cbegin(iterable),
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
constexpr auto sequence(dindex iterable) {

    struct iterator {
        dindex i;

        bool operator!=(const iterator &rhs) const { return i != rhs.i; }

        /** Increase index and iterator at once */
        void operator++() { ++i; }

        /** Tie them together for returning */
        auto operator*() const { return i; }
    };
    struct iterable_wrapper {
        dindex iterable;
        auto begin() { return iterator{iterable}; }
        auto end() { return iterator{iterable + 1}; }
    };
    return iterable_wrapper{std::forward<dindex>(iterable)};
}

/** Helper method to run over a range
 *
 * Usage:
 * for (auto i : sequence(r)) {}
 *
 * with r a range that can be accessed with r[0] and r[1]
 *
 * @note sequence({2,4}) will produce { 2, 3, 4 }
 *
 **/
template <
    typename array_type,
    typename = std::enable_if_t<
        std::conditional_t<std::is_array_v<array_type>, std::extent<array_type>,
                           std::tuple_size<array_type>>::value == 2U>>
constexpr auto sequence(array_type iterable) {

    struct iterator {
        size_t i;
        size_t end;

        bool operator!=(const iterator &rhs) const { return i != rhs.end; }

        /** Increase index and iterator at once */
        void operator++() { ++i; }

        /** Tie them together for returning */
        auto operator*() const { return i; }
    };
    struct iterable_wrapper {
        array_type _iterable;
        auto begin() { return iterator{_iterable[0], _iterable[1]}; }
        auto end() { return iterator{_iterable[1] + 1, _iterable[1] + 1}; }
    };
    return iterable_wrapper{std::forward<array_type>(iterable)};
}

}  // namespace detray
