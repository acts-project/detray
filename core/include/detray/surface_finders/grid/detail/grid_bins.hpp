/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::bins {

/// @brief Bin with a single entry
template <typename entry_t>
class single : public detray::ranges::single_view<entry_t> {

    using base_type = detray::ranges::single_view<entry_t>;

    public:
    using entry_type = entry_t;
    using base_type::base_type;
    using base_type::operator=;

    /// Default constructor initializer the bin with an invalid value
    DETRAY_HOST_DEVICE constexpr single()
        : base_type{detail::invalid_value<entry_t>()} {}

    /// @returns the storage capacity of this bin
    DETRAY_HOST_DEVICE
    constexpr dindex capacity() noexcept { return 1u; }

    /// Add a new entry to the bin
    template <typename E = entry_t>
    DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) noexcept {
        *(*this) = std::forward<E>(entry);
    }

    /// @returns an initialized bin in the backend storage
    DETRAY_HOST_DEVICE
    constexpr auto init(entry_t entry = detail::invalid_value<entry_t>())
        -> single& {
        *(*this) = entry;
        return *this;
    }
};

/// @brief Bin that holds a collection of entries.
///
/// Keeps an additional counter to track the number of entries per bin.
/// The bin capaciy is static.
template <typename entry_t, std::size_t N>
class static_array
    : public detray::ranges::view_interface<static_array<entry_t, N>> {

    using bin_view_t = detray::ranges::subrange<std::array<entry_t, N>>;
    using const_bin_view_t =
        detray::ranges::subrange<const std::array<entry_t, N>>;

    public:
    using entry_type = entry_t;

    /// Default constructor initializer the bin with an invalid value
    DETRAY_HOST_DEVICE constexpr static_array() { init(); };
    constexpr static_array(const static_array& other) = default;
    constexpr static_array(static_array&& other) = default;
    static_array& operator=(const static_array& other) = default;

    /// @returns view iterator over bin content in start or end position
    /// @{
    DETRAY_HOST_DEVICE
    auto begin() { return detray::ranges::begin(view()); }
    DETRAY_HOST_DEVICE
    auto begin() const { return detray::ranges::cbegin(view()); }
    DETRAY_HOST_DEVICE
    auto end() { return detray::ranges::end(view()); }
    DETRAY_HOST_DEVICE
    auto end() const { return detray::ranges::cend(view()); }
    /// @}

    /// @returns the number of entries in this bin - const
    DETRAY_HOST_DEVICE
    constexpr dindex size() const { return m_size; }

    /// The storage capacity of this bin
    DETRAY_HOST_DEVICE
    constexpr dindex capacity() const noexcept {
        return static_cast<dindex>(m_content.size());
    }

    /// Add a new entry to the bin
    template <typename E = entry_t>
    DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) noexcept {
        m_content[m_size] = std::forward<E>(entry);
        ++m_size;
    }

    /// Initilialize with a single entry filled
    ///
    /// @returns Access to the initialized bin
    DETRAY_HOST_DEVICE
    constexpr auto init(entry_t entry = detail::invalid_value<entry_t>())
        -> static_array& {

        // Initialize the storage element
        for (auto& e : m_content) {
            e = detail::invalid_value<entry_t>();
        }

        if (entry == detail::invalid_value<entry_t>()) {
            m_size = 0u;
            return *this;
        }

        // The bin has at least a capacity of 1
        m_content[0] = entry;
        m_size = 1u;

        return *this;
    }

    /// Initilialize from an entire bin content @param content.
    ///
    /// @returns Access to the initialized bin
    DETRAY_HOST_DEVICE
    constexpr auto init(std::array<entry_t, N> content) -> static_array& {

        m_content = content;
        m_size = 0u;
        for (const auto& entry : m_content) {
            if (entry != detail::invalid_value<entry_t>()) {
                ++m_size;
            }
        }

        return *this;
    }

    private:
    /// @returns the subrange on the valid bin content - const
    DETRAY_HOST_DEVICE constexpr auto view() const -> const_bin_view_t {
        return const_bin_view_t{m_content, dindex_range{0u, m_size}};
    }

    /// @returns the subrange on the valid bin content
    DETRAY_HOST_DEVICE constexpr auto view() -> bin_view_t {
        return bin_view_t{m_content, dindex_range{0u, m_size}};
    }

    /// Number of valid elements in the bin
    dindex m_size{0u};
    /// Bin entry container
    std::array<entry_t, N> m_content{};
};

/// @brief Bin that views a collection of entries it does not own.
///
/// Used if the bin capacity is not static.
template <typename container_t>
class dynamic_array
    : public detray::ranges::view_interface<dynamic_array<container_t>> {

    using base_type =
        detray::ranges::view_interface<dynamic_array<container_t>>;
    using bin_view_t = detray::ranges::subrange<container_t>;
    using const_bin_view_t = detray::ranges::subrange<const container_t>;

    public:
    using entry_type = typename container_t::value_type;

    /// Default constructor initializer the bin with an invalid value
    DETRAY_HOST_DEVICE constexpr dynamic_array() { init(); };

    /// Construct from an externally owned container of bin content
    /// @param bin_storage with an explict @param capacity and @param size , as
    /// well as @param offset into the container
    dynamic_array(container_t& bin_storage, const dindex capacity,
                  const dindex offset = 0u, const dindex size = 0u)
        : m_global_storage{&bin_storage},
          m_size{size},
          m_capacity{capacity},
          m_offset{offset} {}

    /// @returns view iterator over bin content in start or end position
    /// @{
    DETRAY_HOST_DEVICE auto begin() { return detray::ranges::begin(view()); }
    DETRAY_HOST_DEVICE auto begin() const {
        return detray::ranges::cbegin(view());
    }
    DETRAY_HOST_DEVICE auto end() { return detray::ranges::end(view()); }
    DETRAY_HOST_DEVICE auto end() const { return detray::ranges::cend(view()); }
    /// @}

    /// @returns the number of entries in this bin - const
    DETRAY_HOST_DEVICE
    constexpr dindex size() const { return m_size; }

    /// The storage capacity of this bin
    DETRAY_HOST_DEVICE
    constexpr dindex capacity() const noexcept { return m_capacity; }

    /// Add a new entry to the bin
    /// @note This does not check the state of the containter it points to!!!
    template <typename E = entry_type>
    DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) noexcept {
        (*m_global_storage)[m_offset + m_size] = std::forward<E>(entry);
        ++m_size;
    }

    /// @note The bin capacity has to be set correctly before calling this
    /// method
    /// @returns Access to an initialized bin in the backend storage
    DETRAY_HOST_DEVICE
    constexpr auto init(entry_type entry = detail::invalid_value<entry_type>())
        -> dynamic_array& {
        if (capacity() == 0u) {
            return *this;
        }

        // Initialize the storage element
        std::fill(this->begin(), this->end(),
                  detail::invalid_value<entry_type>());

        if (entry == detail::invalid_value<entry_type>()) {
            m_size = 0u;
            return *this;
        }

        // The bin has at least a capacity of 1
        this->front() = entry;
        m_size = 1u;

        return *this;
    }

    /// Initilialize from an entire bin content @param content.
    ///
    /// @returns Access to the initialized bin
    DETRAY_HOST_DEVICE
    constexpr auto init(const container_t& content) -> dynamic_array& {

        m_size = 0u;
        for (dindex i{m_offset}; i < m_offset + m_capacity; ++i) {
            if (content[i] != detail::invalid_value<entry_type>()) {
                m_global_storage[m_offset + i] = content[i];
                ++m_size;
            }
        }

        return *this;
    }

    private:
    /// @returns the subrange on the valid bin content - const
    DETRAY_HOST_DEVICE auto view() const -> const_bin_view_t {
        return const_bin_view_t{*m_global_storage,
                                dindex_range{m_offset, m_offset + m_size}};
    }

    /// @returns the subrange on the valid bin content
    DETRAY_HOST_DEVICE auto view() -> bin_view_t {
        return bin_view_t{*m_global_storage,
                          dindex_range{m_offset, m_offset + m_size}};
    }

    /// Pointer to the global bin storage that is not owned by this class
    container_t* m_global_storage{nullptr};
    /// Current number, capacity and offset of bin entries
    dindex m_size{0u}, m_capacity{0u}, m_offset{0u};
};
/// @}

}  // namespace detray::bins
