/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <cassert>
#include <memory>
#include <type_traits>
#include <utility>

namespace detray::ranges {

/// @brief Range adaptor that indexes the elements of a range using another
///        index range.
///
/// @tparam range_itr_t the iterator type of the enumerated range
/// @tparam sequence_itr_t range of indices.
///
/// @note Does not take ownership of the range it operates on. Its lifetime
/// needs to be guranteed throughout iteration or between iterations with the
/// same enumerate instance.
/// @note Is not fit for lazy evaluation.
template <typename range_itr_t, typename sequence_itr_t>
class pick_view : public detray::ranges::view_interface<
                      pick_view<range_itr_t, sequence_itr_t>> {

    private:
    using index_t = typename std::iterator_traits<sequence_itr_t>::value_type;
    using value_t = typename std::iterator_traits<range_itr_t>::value_type;
    using difference_t =
        typename std::iterator_traits<range_itr_t>::difference_type;

    /// @brief Nested iterator to randomly index the elements of a range.
    ///
    /// The indices by which to reference the range are obtained by a dedicated
    /// index range
    ///
    /// @todo Add Comparability to fulfill random access iterator traits once
    ///       needed.
    struct iterator {

        using difference_type =
            typename std::iterator_traits<range_itr_t>::difference_type;
        using value_type = value_t;
        using pointer = typename std::iterator_traits<range_itr_t>::pointer;
        using reference = typename std::iterator_traits<range_itr_t>::reference;
        using iterator_category =
            typename std::iterator_traits<sequence_itr_t>::iterator_category;

        static_assert(std::is_convertible_v<index_t, difference_type>,
                      "Given sequence cannot be "
                      "used to index elements of range.");

        /// @returns true if we reach end of sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator==(const iterator &rhs) const -> bool {
            return (m_seq_iter == rhs.m_seq_iter);
        }

        /// Determine whether we reach end of range - const
        DETRAY_HOST_DEVICE constexpr auto operator!=(const iterator &rhs) const
            -> bool {
            return (m_seq_iter != rhs.m_seq_iter);
        }

        /// Increment iterator and index in lockstep
        DETRAY_HOST_DEVICE constexpr auto operator++() -> iterator & {
            ++m_seq_iter;
            if (m_seq_iter != m_seq_end) {
                m_range_iter =
                    m_range_begin + static_cast<difference_type>(*(m_seq_iter));
            }
            return *this;
        }

        /// Decrement iterator and index in lockstep
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::bidirectional_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator--() -> iterator & {
            m_range_iter =
                m_range_begin + static_cast<difference_type>(*(--m_seq_iter));
            return *this;
        }

        /// @returns iterator and index together
        DETRAY_HOST_DEVICE auto operator*() {
            return std::pair<
                const typename std::iterator_traits<sequence_itr_t>::value_type,
                const typename std::iterator_traits<range_itr_t>::value_type &>(
                *m_seq_iter, *m_range_iter);
        }

        /// @returns an iterator and index position advanced by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator+(
            const difference_type j) const -> iterator {
            auto seq_iter = detray::ranges::next(m_seq_iter, j);
            return {m_range_begin + static_cast<difference_type>(*seq_iter),
                    m_range_begin, seq_iter, m_seq_end};
        }

        /// @returns an iterator and index position advanced by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-(
            const difference_type j) const -> iterator {
            return *this + -j;
        }

        /// @returns the positional difference between two iterations
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-(const iterator &other) const
            -> difference_type {
            return m_seq_iter - other.m_seq_iter;
        }

        /// @returns advance this iterator state by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
            -> iterator & {
            detray::ranges::advance(m_seq_iter, j);
            m_range_begin += static_cast<difference_type>(*m_seq_iter);
            return *this;
        }

        /// @returns advance this iterator state by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-=(const difference_type j)
            -> iterator & {
            return *this += -j;
        }
        range_itr_t m_range_iter, m_range_begin;
        sequence_itr_t m_seq_iter, m_seq_end;
    };

    range_itr_t m_range_begin, m_range_end;
    sequence_itr_t m_seq_begin, m_seq_end;

    public:
    using iterator_t = iterator;

    /// Default constructor
    pick_view() = default;

    /// Construct from a @param range that will be enumerated beginning at 0
    template <
        typename range_t, typename sequence_t,
        std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true,
        std::enable_if_t<detray::ranges::range_v<sequence_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr pick_view(range_t &&range, sequence_t &&seq)
        : m_range_begin{detray::ranges::begin(std::forward<range_t>(range))},
          m_range_end{detray::ranges::end(std::forward<range_t>(range))},
          m_seq_begin{detray::ranges::cbegin(std::forward<sequence_t>(seq))},
          m_seq_end{detray::ranges::cend(std::forward<sequence_t>(seq))} {}

    /// Copy constructor
    DETRAY_HOST_DEVICE
    constexpr pick_view(const pick_view &other)
        : m_range_begin{other.m_range_begin},
          m_range_end{other.m_range_end},
          m_seq_begin{other.m_seq_begin},
          m_seq_end{other.m_seq_end} {}

    /// Default destructor
    DETRAY_HOST_DEVICE ~pick_view() {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    pick_view &operator=(const pick_view &other) {
        m_range_begin = other.m_range_begin;
        m_range_end = other.m_range_end;
        m_seq_begin = other.m_seq_begin;
        m_seq_end = other.m_seq_end;
        return *this;
    }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator {
        return {m_range_begin + static_cast<difference_t>(*(m_seq_begin)),
                m_range_begin, m_seq_begin, m_seq_end};
    }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator {
        return {m_range_begin + static_cast<difference_t>(*(m_seq_begin)),
                m_range_begin, m_seq_begin, m_seq_end};
    }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() -> iterator {
        return {m_range_begin, m_range_begin, m_seq_end, m_seq_end};
    }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const -> const typename iterator_t::value_type * {
        return std::addressof(
            *(m_range_begin + static_cast<difference_t>(*m_seq_begin)));
    }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() -> typename iterator_t::value_type * {
        return std::addressof(
            *(m_range_begin + static_cast<difference_t>(*m_seq_begin)));
    }

    /// @returns the number of elements in the underlying range. Simplified
    /// implementation compared to @c view_interface.
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept {
        return detray::ranges::distance(m_seq_begin, m_seq_end);
    }

    /// @returns access to the first element value in the range, including the
    /// corresponding index.
    DETRAY_HOST_DEVICE
    constexpr auto front() noexcept {
        return std::pair<index_t, value_t &>(
            *(m_seq_begin),
            *detray::ranges::next(m_range_begin,
                                  static_cast<difference_t>(*m_seq_begin)));
    }

    /// @returns access to the last element value in the range, including the
    /// corresponding index.
    DETRAY_HOST_DEVICE
    constexpr auto back() noexcept {
        index_t last_idx{*detray::ranges::next(
            m_seq_begin, static_cast<difference_t>((size() - 1u)))};
        return std::pair<index_t, value_t &>(
            last_idx, *detray::ranges::next(
                          m_range_begin, static_cast<difference_t>(last_idx)));
    }

    DETRAY_HOST_DEVICE
    constexpr auto operator[](const dindex i) const {
        index_t last_idx{
            *detray::ranges::next(m_seq_begin, static_cast<difference_t>(i))};
        return std::pair<index_t, value_t &>(
            last_idx, *detray::ranges::next(
                          m_range_begin, static_cast<difference_t>(last_idx)));
    }
};

namespace views {

template <typename range_itr_t, typename sequence_itr_t>
struct pick : public pick_view<range_itr_t, sequence_itr_t> {

    using base_type = pick_view<range_itr_t, sequence_itr_t>;

    pick() = default;

    template <
        typename range_t, typename sequence_t,
        std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true,
        std::enable_if_t<detray::ranges::range_v<sequence_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr pick(range_t &&range, sequence_t &&seq)
        : base_type(std::forward<range_t>(range),
                    std::forward<sequence_t>(seq)) {}
};

// deduction guides

template <typename range_t, typename sequence_t,
          std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true,
          std::enable_if_t<detray::ranges::range_v<sequence_t>, bool> = true>
pick(range_t &&range, sequence_t &&seq)
    -> pick<detray::ranges::iterator_t<range_t>,
            detray::ranges::const_iterator_t<sequence_t>>;

}  // namespace views

}  // namespace detray::ranges
