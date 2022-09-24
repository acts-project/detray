/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/ranges/subrange.hpp"

// System include(s)
#include <tuple>
#include <type_traits>

namespace detray::ranges {

/// @brief Enumerates the elements of a range on the fly.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_itr_t the iterator type of the enumerated range
/// @tparam incr_t a type that can be incremented in lockstep with the
///         iterator 'range_itr_t'.
template <typename range_itr_t, typename sequence_itr_t>
class pick_view : public detray::ranges::view_interface<
                      pick_view<range_itr_t, sequence_itr_t>> {

    private:
    /// @brief Nested iterator to enumerate the elements of a range.
    ///
    /// The enumeration is done by incrementing an index in lockstep with a
    /// wrapped iterator of a range. Index and current iterator value are
    /// returned using structured binding.
    struct iterator {

        using difference_type =
            typename std::iterator_traits<range_itr_t>::difference_type;
        using value_type =
            typename std::iterator_traits<range_itr_t>::value_type;
        using pointer = typename std::iterator_traits<range_itr_t>::pointer;
        using reference = typename std::iterator_traits<range_itr_t>::reference;
        using iterator_category =
            typename std::iterator_traits<sequence_itr_t>::iterator_category;

        /// Determine whether we reach end of range
        DETRAY_HOST_DEVICE
        constexpr auto operator!=(const iterator &rhs) const -> bool {
            return m_seq_iter != rhs.m_seq_iter;
        }

        /// Increment iterator and index in lockstep
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            m_range_iter = m_range_begin + *(++m_seq_iter);
            return *this;
        }

        /// @return iterator and index together
        DETRAY_HOST_DEVICE
        auto operator*() { return std::tie(*m_seq_iter, *m_range_iter); }

        /// Advance the iterator and index position by @param j.
        DETRAY_HOST_DEVICE
        constexpr auto operator+(const difference_type j) const -> iterator {
            return {std::next(m_range_begin, *(m_seq_iter + j)), m_range_begin,
                    m_seq_iter + j, m_seq_end};
        }

        range_itr_t m_range_iter, m_range_begin;
        sequence_itr_t m_seq_iter, m_seq_end;
    };

    range_itr_t m_range_begin, m_range_end;
    sequence_itr_t m_seq_begin, m_seq_end;

    public:
    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    pick_view() = default;

    /// Construct from a @param range that will be enumerated beginning at 0
    template <typename range_t, typename sequence_t>
    DETRAY_HOST_DEVICE constexpr pick_view(range_t &&range, sequence_t &&seq)
        : m_range_begin{detray::ranges::begin(std::forward<range_t>(range))},
          m_range_end{detray::ranges::end(std::forward<range_t>(range))},
          m_seq_begin{detray::ranges::begin(std::forward<sequence_t>(seq))},
          m_seq_end{detray::ranges::end(std::forward<sequence_t>(seq))} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator {
        return {m_range_begin + *(m_seq_begin), m_range_begin, m_seq_begin,
                m_seq_end};
    }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator {
        return {m_range_begin + *(m_seq_begin), m_range_begin, m_seq_begin,
                m_seq_end};
    }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() -> iterator {
        return {m_range_end, m_range_begin, m_seq_end, m_seq_end};
    }
};

namespace views {

template <typename range_itr_t, typename sequence_itr_t>
struct pick : public pick_view<range_itr_t, sequence_itr_t> {

    using base_type = pick_view<range_itr_t, sequence_itr_t>;

    pick() = default;

    template <typename range_t, typename sequence_t>
    DETRAY_HOST_DEVICE constexpr pick(range_t &&range, sequence_t &&seq)
        : base_type(std::forward<range_t>(range),
                    std::forward<sequence_t>(seq)) {}
};

// deduction guides

template <typename range_t, typename sequence_t>
pick(range_t &&range, sequence_t &&seq)
    -> pick<detray::ranges::iterator_t<range_t>,
            detray::ranges::iterator_t<sequence_t>>;

}  // namespace views

}  // namespace detray::ranges