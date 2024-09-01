/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/ranges/subrange.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray::ranges {

/// @brief Enumerates the elements of a range on the fly.
///
/// @tparam range_itr_t the iterator type of the enumerated range
/// @tparam incr_t a type that can be incremented in lockstep with the
///         iterator 'range_itr_t'.
///
/// @note Does not take ownership of the range it operates on. Its lifetime
/// needs to be guranteed throughout iteration or between iterations with the
/// same enumerate instance.
/// @note Is not fit for lazy evaluation.
template <typename range_itr_t, typename incr_t = std::size_t>
class enumerate_view : public detray::ranges::view_interface<
                           enumerate_view<range_itr_t, incr_t>> {

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
            typename std::iterator_traits<range_itr_t>::iterator_category;

        constexpr iterator() = default;

        DETRAY_HOST_DEVICE
        iterator(range_itr_t iter, incr_t offset = 0)
            : m_iter(iter), m_i{offset}, m_offset{offset} {};

        /// @returns true if we reach end of sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator==(const iterator &rhs) const -> bool {
            return (m_iter == rhs.m_iter);
        }

        /// @returns true if the wrapped iterators are not the same.
        DETRAY_HOST_DEVICE constexpr auto operator!=(const iterator &rhs) const
            -> bool {
            return (m_iter != rhs.m_iter);
        }

        /// Increment iterator and index in lockstep
        DETRAY_HOST_DEVICE constexpr auto operator++() -> iterator & {
            ++m_i;
            ++m_iter;
            return *this;
        }

        /// Increment iterator and index in lockstep
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::bidirectional_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator--() -> iterator & {
            --m_i;
            --m_iter;
            return *this;
        }

        /// @returns iterator and index together
        DETRAY_HOST_DEVICE constexpr auto operator*() const {
            return std::pair<incr_t, const value_type &>(m_i, *m_iter);
        }

        /// @returns an iterator and index position advanced by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator+(
            const difference_type j) const -> iterator {
            return {m_iter + j, m_i + static_cast<incr_t>(j)};
        }

        /// @returns an iterator and index position advanced by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-(
            const difference_type j) const -> iterator {
            return {m_iter - j, m_i - static_cast<incr_t>(j)};
        }

        /// @returns the positional difference between two iterators
        /// (independent from their enumeration of the range values)
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-(const iterator &other) const
            -> difference_type {
            return m_iter - other.m_iter;
        }

        /// @returns advance this iterator state by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
            -> iterator & {
            m_iter += j;
            m_i += static_cast<incr_t>(j);
            return *this;
        }

        /// @returns advance this iterator state by @param j.
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator-=(const difference_type j)
            -> iterator & {
            m_iter -= j;
            m_i -= static_cast<incr_t>(j);
            return *this;
        }

        /// @returns the value and index at a given position - const
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator[](
            const difference_type i) const {
            // check narrowing
            const incr_t index{m_offset + static_cast<incr_t>(i)};
            return std::pair<incr_t, const value_type &>(index, m_iter[i]);
        }

        /// @returns the value and index at a given position - const
        template <typename I = range_itr_t,
                  std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                                   bool> = true>
        DETRAY_HOST_DEVICE constexpr auto operator[](const difference_type i) {
            // check narrowing
            const incr_t index{m_offset + static_cast<incr_t>(i)};
            return std::pair<incr_t, const value_type &>(index, m_iter[i]);
        }

        range_itr_t m_iter{};
        incr_t m_i{0}, m_offset{0};
    };

    iterator m_begin, m_end;

    public:
    using iterator_t = iterator;

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    enumerate_view() = default;

    /// Construct from a @param range that will be enumerated beginning at 0
    template <typename range_t,
              std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(range_t &&rng)
        : m_begin{detray::ranges::begin(std::forward<range_t>(rng))},
          m_end{detray::ranges::end(std::forward<range_t>(rng)),
                static_cast<dindex>(rng.size())} {}

    /// Construct from a @param range that will be enumerated beginning at
    /// @param start.
    template <typename range_t,
              std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr enumerate_view(range_t &&rng, dindex start)
        : m_begin{detray::ranges::begin(std::forward<range_t>(rng)), start},
          m_end{detray::ranges::end(std::forward<range_t>(rng)),
                start + static_cast<dindex>(rng.size())} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator { return m_begin; }

    /// @return start position of range on container - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator { return m_begin; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator { return m_end; }
};

namespace views {

template <
    typename range_itr_t, typename incr_t = dindex,
    std::enable_if_t<std::is_convertible_v<typename std::iterator_traits<
                                               range_itr_t>::difference_type,
                                           incr_t>,
                     bool> = true>
struct enumerate : public enumerate_view<range_itr_t, incr_t> {

    using base_type = enumerate_view<range_itr_t, incr_t>;

    constexpr enumerate() = default;

    template <typename range_t,
              std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate(range_t &&rng)
        : base_type(std::forward<range_t>(rng)) {}

    template <typename range_t,
              std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr enumerate(range_t &&rng, incr_t start)
        : base_type(std::forward<range_t>(rng), start) {}

    /// Construct from a @param range and an index range provided by a volume
    /// @param vol.
    template <typename deduced_range_t, typename volume_t,
              typename = typename std::remove_reference_t<volume_t>::volume_def>
    DETRAY_HOST_DEVICE enumerate(deduced_range_t &&range, const volume_t &vol)
        : enumerate(detray::ranges::subrange(
                        std::forward<deduced_range_t>(range), vol),
                    detray::detail::get<0>(vol.full_range())) {}
};

// deduction guides

template <typename range_t,
          std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t, typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE enumerate(range_t &&range, const volume_t &vol)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t,
          std::enable_if_t<detray::ranges::range_v<range_t>, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng, dindex start)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

}  // namespace views

}  // namespace detray::ranges
