/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// System include(s)
#include <iostream>
#include <tuple>  // std::tie
#include <type_traits>

// Project include(s)
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// @brief Struct that implements a range by providing start and end iterators.
///
/// @tparam container_t the container on which to perform the iteration
template <typename container_t,
          typename = typename std::remove_reference_t<container_t>::value_type>
struct iterator_range {
    using container_iter_type =
        decltype(std::begin(std::declval<container_t>()));
    using value_type = typename container_t::value_type;

    /// Delete default constructor
    iterator_range() = delete;

    /// Always construct from a container and a range
    ///
    /// @param iterable container to iterate over
    /// @param range start and end position for iteration. Must contain two
    /// arithmetic type values
    template <typename range_t>
    DETRAY_HOST_DEVICE iterator_range(const container_t &iterable,
                                      const range_t &range)
        : _start(std::begin(iterable) + detail::get<0>(range)),
          _end(std::begin(iterable) + detail::get<1>(range)) {
        // Observe end of the iterable
        if (std::distance(_end, std::end(iterable)) <= 0) {
            _end = std::end(iterable);
        }
    }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    inline auto begin() -> container_iter_type & { return _start; }

    /// @return end position of range on container.
    DETRAY_HOST_DEVICE
    inline auto end() -> container_iter_type & { return _end; }

    /// Does this describe the same range?
    DETRAY_HOST_DEVICE
    inline auto operator!=(const iterator_range &rhs) const -> bool {
        return _start != rhs._start or _end != rhs._end;
    }

    /// @return element at position i, relative to iterator range.
    DETRAY_HOST_DEVICE
    inline auto operator[](const dindex i) -> value_type & {
        return *(_start + i);
    }

    /// @return element at position i, relative to iterator range - const
    DETRAY_HOST_DEVICE
    inline auto operator[](const dindex i) const -> const value_type & {
        return *(_start + i);
    }

    /// @return the offset of the range start into the container.
    DETRAY_HOST_DEVICE
    inline auto offset(const container_t &iterable) const -> std::size_t {
        return std::distance(_start, std::begin(iterable));
    }

    /// Start and end position of a range
    container_iter_type _start, _end;
};

/// Get an interator range for the constituents of a volume on their container.
///
/// @tparam container_t the container on which to perform the iteration
/// @tparam volume_t the type of volume from which to get the range in the
///         container
///
/// @param iterable reference to the container
/// @param volume the volume
///
/// @returns an iterator range on the container, according to the volumes range
template <typename container_t,
          typename = typename std::remove_reference_t<container_t>::value_type,
          typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_t &iterable,
                                               volume_t &&volume) {
    return iterator_range(
        iterable, volume.template range<typename container_t::value_type>());
}

/// Overload of the range-function for dindex_range
template <typename container_t>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_t &iterable,
                                               const dindex_range &range) {

    return iterator_range(iterable, range);
}

/// Overload of the range-function for a single index
template <typename container_t>
DETRAY_HOST_DEVICE inline constexpr auto range(const container_t &iterable,
                                               const dindex &i) {

    return iterator_range(iterable, dindex_range{i, i + 1});
}

/// @brief Struct that increments a given iterator in lockstep with the
///        iteration index for convenience.
///
/// @tparam container_t the container on which to perform the iteration
/// @tparam volume_t the type of volume from which to get the range in the
///         container
template <typename container_t,
          typename = typename std::remove_reference_t<container_t>::value_type,
          typename = decltype(std::end(std::declval<container_t>()))>
struct enumerator {

    using container_iter_type =
        decltype(std::begin(std::declval<container_t>()));
    using value_type =
        typename std::remove_reference_t<container_t>::value_type;

    /// Start position of iterator
    std::size_t i;
    /// Iterator on container
    container_iter_type iter;

    /// Build from start index and corresponding iterator - rvalue
    DETRAY_HOST_DEVICE
    enumerator(std::size_t start, const container_iter_type &&iterator)
        : i(start), iter(iterator) {}

    /// Build from start index and corresponding iterator - lvalue
    DETRAY_HOST_DEVICE
    enumerator(std::size_t start, const container_iter_type &iterator)
        : i(start), iter(iterator) {}

    /// Determine end of iteration
    DETRAY_HOST_DEVICE
    auto operator!=(const enumerator &rhs) const -> bool {
        return iter != rhs.iter;
    }

    /// Increase index and iterator at once
    DETRAY_HOST_DEVICE
    auto operator++() -> void {
        ++i;
        ++iter;
    }

    /// Tie them together for returning
    DETRAY_HOST_DEVICE
    auto operator*() const { return std::tie(i, *iter); }
};

/// Helper utility to allow indexed enumeration with structured binding
///
/// Usage:
///
/// for (auto [ i, value ] = enumerate(container) ) { ... };
///
/// with 'container' any stl-like container
template <typename container_t,
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

/// Helper utility to allow indexed enumeration with structured binding for
///  a given range
///
/// Usage:
///
/// for (auto [ i, value ] = enumerate(container, range) ) { ... };
///
/// with 'container' any stl-like container and r a range type that is
/// compatible with std::get or a type for which an overload to the range()
/// function exists.
template <typename container_t, typename range_t,
          typename = decltype(std::end(std::declval<container_t>())),
          typename = std::enable_if_t<
              std::is_integral_v<std::remove_reference_t<
                  decltype(std::get<0>(std::declval<range_t>()))>> and
              std::is_integral_v<std::remove_reference_t<
                  decltype(std::get<1>(std::declval<range_t>()))>>>>
DETRAY_HOST_DEVICE constexpr inline auto enumerate(const container_t &iterable,
                                                   range_t &&r) {

    struct iterable_wrapper {

        using container_iter_type =
            decltype(std::begin(std::declval<container_t>()));

        const container_iter_type iter;
        iterator_range<container_t> range_iter;

        DETRAY_HOST_DEVICE
        auto begin() -> enumerator<container_t> {
            return enumerator<container_t>(
                static_cast<std::size_t>(
                    std::distance(iter, range_iter.begin())),
                range_iter.begin());
        }
        DETRAY_HOST_DEVICE
        auto end() -> enumerator<container_t> {
            return enumerator<container_t>(
                static_cast<std::size_t>(
                    std::distance(iter, range_iter.begin())),
                range_iter.end());
        }
    };

    return iterable_wrapper{std::begin(iterable),
                            range(iterable, std::forward<range_t>(r))};
}

template <typename container_t,
          typename = typename std::remove_reference_t<container_t>::value_type,
          typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE constexpr inline auto enumerate(const container_t &iterable,
                                                   volume_t &&vol) {
    return enumerate(iterable, vol.range());
}

/// Helper method to run over an irregular sequence
///
/// Usage:
/// for (auto i : sequence(r)) {}
///
/// with r an index range (integral) that can be accessed with std::get
///
/// @note sequence({2,4}) will produce { 2, 3, 4 }
template <typename container_t>
DETRAY_HOST_DEVICE auto enumerate(const container_t &iterable,
                                  const dvector<dindex> indices) {

    struct iterator {
        using sequence_iter_type =
            decltype(std::begin(std::declval<dvector<dindex>>()));
        using container_iter_type =
            decltype(std::begin(std::declval<container_t>()));
        using value_type =
            typename std::remove_reference_t<container_t>::value_type;

        container_iter_type cont_begin;
        container_iter_type cont_iter;
        sequence_iter_type seq_iter;
        sequence_iter_type seq_end;

        /// Determine whether we reach end of sequence
        DETRAY_HOST_DEVICE
        auto operator!=(const iterator &rhs) const -> bool {
            return (seq_iter != std::next(rhs.seq_end, -1));
        }

        /// Increase index and iterator at once
        DETRAY_HOST_DEVICE
        auto operator++() -> void {
            // Reset container iterator and have it point to the current surface
            // cont_iter = cont_begin;
            // std::cout << "advance: " << *seq_iter << std::endl;
            cont_iter = std::next(cont_begin, *(++seq_iter));

            // std::cout << "new pos: " << std::distance(cont_iter, cont_begin)
            // << std::endl;
            //  Go to the next surface

            // std::cout << "next time: " << *seq_iter << std::endl;

            /*std::cout << "after: " << (*cont_iter).volume() << std::endl;
            std::cout << "advanced by: " << *seq_iter << ", " << *(++seq_iter)
            << ", " << std::abs(static_cast<int>(*seq_iter) -
            static_cast<int>(*(++seq_iter))) << std::endl;*/
        }

        /// Tie them together for returning
        DETRAY_HOST_DEVICE
        auto operator*() { return std::tie(*seq_iter, *cont_iter); }
    };

    /// Wrap up for iteration
    struct iterable_wrapper {
        const container_t &_iterable;
        const dvector<dindex> _indices;

        DETRAY_HOST_DEVICE
        auto begin() {
            // std::cout << "New neighborhood" << std::endl;
            // for (const dindex i : _indices)
            //     std::cout << i << std::endl;

            // std::cout << "size: " << _indices.size() << std::endl;
            return iterator{std::begin(_iterable), std::begin(_iterable),
                            std::begin(_indices), std::end(_indices)};
        }

        DETRAY_HOST_DEVICE
        auto end() {
            return iterator{std::end(_iterable), std::end(_iterable),
                            std::end(_indices), std::end(_indices)};
        }
    };

    return iterable_wrapper{iterable, indices};
}

/// Helper method to (fake a) run over a single entry
///
/// Usage:
/// for (auto i : sequence(j)) {}
///
/// with j an unsinged dindex type
///
/// @note sequence(2) will produce {2}
DETRAY_HOST_DEVICE
constexpr auto sequence(dindex iterable) {

    struct iterator {
        /// Start and end of sequence
        dindex i;
        std::size_t end;

        /// Determine whether we reach end of sequence
        DETRAY_HOST_DEVICE
        auto operator!=(const iterator &rhs) const -> bool {
            return i != rhs.end;
        }

        /// Increase index and iterator at once
        DETRAY_HOST_DEVICE
        auto operator++() -> void { ++i; }

        /// @returns the current index in the sequence
        DETRAY_HOST_DEVICE
        auto operator*() const -> dindex { return i; }
    };

    /// Wrap up for iteration
    struct iterable_wrapper {
        dindex iterable;

        DETRAY_HOST_DEVICE
        auto begin() { return iterator{iterable, iterable + 1}; }

        DETRAY_HOST_DEVICE
        auto end() { return iterator{iterable + 1, iterable + 1}; }
    };

    return iterable_wrapper{iterable};
}

/// Helper method to run over a range
///
/// Usage:
/// for (auto i : sequence(r)) {}
///
/// with r an index range (integral) that can be accessed with std::get
///
/// @note sequence({2,4}) will produce { 2, 3, 4 }
template <
    typename range_t,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<0>(std::declval<range_t>()))>>>,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<1>(std::declval<range_t>()))>>>>
DETRAY_HOST_DEVICE constexpr auto sequence(range_t &&range) {

    struct iterator {
        /// Start and end of sequence
        std::size_t i;
        std::size_t end;

        /// Determine whether we reach end of sequence
        DETRAY_HOST_DEVICE
        auto operator!=(const iterator &rhs) const -> bool {
            return i < rhs.end;
        }

        /// Increase index and iterator at once
        DETRAY_HOST_DEVICE
        auto operator++() -> void { ++i; }

        /// @returns the current index in the sequence
        DETRAY_HOST_DEVICE
        auto operator*() const -> std::size_t { return i; }
    };

    /// Wrap up for iteration
    struct iterable_wrapper {
        range_t _range;

        DETRAY_HOST_DEVICE
        auto begin() {
            return iterator{std::get<0>(_range), std::get<1>(_range)};
        }
        DETRAY_HOST_DEVICE
        auto end() {
            return iterator{std::get<1>(_range) + 1, std::get<1>(_range) + 1};
        }
    };

    return iterable_wrapper{std::forward<range_t>(range)};
}

}  // namespace detray
