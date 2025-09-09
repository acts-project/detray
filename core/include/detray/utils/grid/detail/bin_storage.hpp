/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#ifndef DETRAY_COMPILE_VITIS
#include "detray/core/detail/container_buffers.hpp"
#endif
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/grid/detail/grid_bins.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

/// @brief bin data state of a grid
///
/// Can be data-owning or not. Does not contain the data of the axes,
/// as that is managed by the multi-axis type directly.
template <bool is_owning, typename bin_t, typename containers>
class bin_storage : public detray::ranges::view_interface<
                        bin_storage<is_owning, bin_t, containers>> {

    template <typename T>
    using vector_t = typename containers::template vector_type<T>;
    using bin_range_t =
        std::conditional_t<is_owning, vector_t<bin_t>,
                           detray::ranges::subrange<vector_t<bin_t>>>;

    public:
    /// Bin type: single or static_array
    using bin_type = bin_t;
    /// Backend storage type for the grid
    using bin_container_type = vector_t<bin_t>;

    // Vecmem based view type
    using view_type = dvector_view<bin_type>;
    using const_view_type = dvector_view<const bin_type>;

    // Vecmem based buffer type
#ifndef DETRAY_COMPILE_VITIS
    using buffer_type = dvector_buffer<bin_type>;
#endif

    /// Default constructor
    bin_storage() = default;
    /// Copy constructor
    bin_storage(const bin_storage&) = default;
    /// Move constructor
    bin_storage(bin_storage&&) = default;

    /// Copy assignment
    bin_storage& operator=(const bin_storage&) = default;
    /// Move assignment
    bin_storage& operator=(bin_storage&&) = default;

    /// Construct containers using a memory resources
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST bin_storage(vecmem::memory_resource& resource)
        : m_bin_data(&resource) {}
#endif // DETRAY_COMPILE_VITIS

    /// Construct grid data from containers - move
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type&& bin_data)
        : m_bin_data(std::move(bin_data)) {}

    /// Construct the non-owning type from the @param offset into the global
    /// container @param bin_data and the number of bins @param size
    template <bool owner = is_owning, std::enable_if_t<!owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type& bin_data, dindex offset,
                                   dindex size)
        : m_bin_data(bin_data, dindex_range{offset, offset + size}) {}

    /// Construct the non-owning type from the @param offset into the global
    /// container @param bin_data and the number of bins @param size
    template <bool owner = is_owning, std::enable_if_t<!owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(const bin_container_type& bin_data,
                                   dindex offset, dindex size)
        : m_bin_data(bin_data, dindex_range{offset, offset + size}) {}

    /// Construct bin storage from its vecmem view
    template <typename view_t,
              typename std::enable_if_t<detail::is_device_view<view_t>::value ,
                                        bool> = true>
    DETRAY_HOST_DEVICE bin_storage(const view_t& view) : m_bin_data(view) {}

    /// begin and end of the bin range
    /// @{
    DETRAY_HOST_DEVICE
    auto begin() { return detray::ranges::begin(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto begin() const { return detray::ranges::cbegin(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto end() { return detray::ranges::end(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto end() const { return detray::ranges::cend(m_bin_data); }
    /// @}

    /// @returns the vecmem view of the bin storage
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return detray::get_data(m_bin_data);
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns the vecmem view of the bin storage - const
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return detray::get_data(m_bin_data);
    }
#endif // DETRAY_COMPILE_VITIS

    private:
    /// Container that holds all bin data when owning or a view into an
    /// externally owned container
    bin_range_t m_bin_data{};
};

/// Facade/wrapper for the data containers of the dynamic bin storage to fit in
/// the grid collection
template <typename bin_t, typename containers>
struct dynamic_bin_container {

    template <typename T>
    using vector_t = typename containers::template vector_type<T>;
    using bin_data_t = typename bin_t::data;

    vector_t<bin_data_t> bins{};
    vector_t<typename bin_t::entry_type> entries{};

    // Vecmem based view type
    using view_type = dmulti_view<dvector_view<bin_data_t>,
                                  dvector_view<typename bin_t::entry_type>>;
    using const_view_type =
        dmulti_view<dvector_view<const bin_data_t>,
                    dvector_view<const typename bin_t::entry_type>>;

    // Vecmem based buffer type
#ifndef DETRAY_COMPILE_VITIS
    using buffer_type =
        dmulti_buffer<dvector_buffer<bin_data_t>,
                      dvector_buffer<typename bin_t::entry_type>>;
#endif

    constexpr dynamic_bin_container() = default;
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    dynamic_bin_container(vecmem::memory_resource* resource)
        : bins{resource}, entries{resource} {}
#endif // DETRAY_COMPILE_VITIS
    dynamic_bin_container(const dynamic_bin_container& other) = default;
    dynamic_bin_container(dynamic_bin_container&& other) = default;

    dynamic_bin_container& operator=(const dynamic_bin_container&) = default;
    dynamic_bin_container& operator=(dynamic_bin_container&&) = default;

    /// Device-side construction from a vecmem based view type
    template <typename view_t,
              typename std::enable_if_t<detail::is_device_view<view_t>::value ,
                                        bool> = true>
    DETRAY_HOST_DEVICE dynamic_bin_container(view_t& view)
        : bins(detail::get<0>(view.m_view)),
          entries(detail::get<1>(view.m_view)) {}

    /// Insert bin data at the end
#ifndef DETRAY_COMPILE_VITIS
    template <typename grid_bin_range_t>
    DETRAY_HOST void append(const grid_bin_range_t& grid_bins) {

        const auto& g_bins = grid_bins.bin_data();
        const auto& g_entries = grid_bins.entry_data();

        bins.insert(bins.end(), g_bins.begin(), g_bins.end());

        // Update the bin offsets
        dindex_range new_bins{
            static_cast<dindex>(bins.size() - grid_bins.size()),
            static_cast<dindex>(bins.size())};
        for (auto& bin : detray::ranges::subrange{bins, new_bins}) {
            bin.update_offset(entries.size());
        }

        entries.insert(entries.end(), g_entries.begin(), g_entries.end());
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns a vecmem view on the bin data - non-const
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST auto get_data() -> view_type {
        return view_type{detray::get_data(bins), detray::get_data(entries)};
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns a vecmem view on the bin data - const
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    auto get_data() const -> const_view_type {
        return const_view_type{detray::get_data(bins),
                               detray::get_data(entries)};
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns the number of bins
    DETRAY_HOST_DEVICE
    std::size_t size() const { return bins.size(); }

    /// Clear out all data
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    void clear() {
        bins.clear();
        entries.clear();
    }
#endif // DETRAY_COMPILE_VITIS
};

/// @brief bin data state of a grid with dynamic bin capacities
///
/// Can be data-owning or not. Does not contain the data of the axes,
/// as that is managed by the multi-axis type directly.
template <bool is_owning, typename entry_t, typename containers>
class bin_storage<is_owning, detray::bins::dynamic_array<entry_t>, containers>
    : public detray::ranges::view_interface<bin_storage<
          is_owning, detray::bins::dynamic_array<entry_t>, containers>> {

    template <typename T>
    using vector_t = typename containers::template vector_type<T>;
    using bin_t = detray::bins::dynamic_array<entry_t>;
    using bin_data_t = typename bin_t::data;
    using bin_range_t =
        std::conditional_t<is_owning, vector_t<bin_data_t>,
                           detray::ranges::subrange<vector_t<bin_data_t>>>;
    using bin_iterator_t = typename detray::ranges::iterator_t<bin_range_t>;
    using const_bin_iterator_t =
        typename detray::ranges::const_iterator_t<bin_range_t>;

    using entry_range_t =
        std::conditional_t<is_owning, vector_t<entry_t>,
                           detray::ranges::subrange<vector_t<entry_t>>>;

    /// Iterator adapter that makes sure the bin storage returns a correctly
    /// initialized bin instance
    template <typename bin_itr_t>
    struct iterator_adapter {
        using difference_type =
            typename std::iterator_traits<bin_itr_t>::difference_type;
        using value_type = bin_t;
        using pointer = bin_t*;
        using reference = bin_t&;
        using iterator_category =
            typename std::iterator_traits<bin_itr_t>::iterator_category;

        using data_ptr_t =
            std::conditional_t<std::is_same<bin_itr_t, const_bin_iterator_t>::value ,
                               const entry_t*, entry_t*>;

        DETRAY_HOST_DEVICE
        iterator_adapter(bin_itr_t&& itr, data_ptr_t entry_data)
            : m_entry_data{entry_data}, m_itr{std::move(itr)} {}

        /// Wrap iterator functionality
        /// @{
        DETRAY_HOST_DEVICE bool operator==(
            const iterator_adapter& other) const {
            return m_itr == other.m_itr;
        }
        DETRAY_HOST_DEVICE bool operator!=(
            const iterator_adapter& other) const {
            return m_itr != other.m_itr;
        }
        DETRAY_HOST_DEVICE iterator_adapter& operator++() {
            ++m_itr;
            return *this;
        }
        DETRAY_HOST_DEVICE iterator_adapter& operator--() {
            --m_itr;
            return *this;
        }
        DETRAY_HOST_DEVICE
        difference_type operator-(const iterator_adapter& other) {
            return m_itr - other.m_itr;
        }
        DETRAY_HOST_DEVICE
        iterator_adapter operator-(difference_type i) {
            return {m_itr - i, m_entry_data};
        }
        DETRAY_HOST_DEVICE
        iterator_adapter operator+(difference_type i) {
            return {m_itr + i, m_entry_data};
        }
        DETRAY_HOST_DEVICE
        constexpr decltype(auto) operator[](const difference_type i) const {
            return *(*this + i);
        }
        DETRAY_HOST_DEVICE
        constexpr decltype(auto) operator[](const difference_type i) {
            return *(*this + i);
        }
        /// @}

        /// Construct and @returns a bin on the fly
        /// @{
#ifndef DETRAY_COMPILE_VITIS
        DETRAY_HOST_DEVICE
        constexpr auto operator*() const {
            return detray::bins::dynamic_array{m_entry_data, *m_itr};
        }
        DETRAY_HOST_DEVICE
        constexpr auto operator*() {
            return detray::bins::dynamic_array{m_entry_data, *m_itr};
        }
        /// @}
#endif // DETRAY_COMPILE_VITIS

        private:
        /// Access to the bin content
        data_ptr_t m_entry_data;
        /// Iterator over bin data from which to construct a bin instance
        bin_itr_t m_itr;
    };

    public:
    /// Bin type: dynamic_array
    using bin_type = bin_t;
    /// Backend storage type for the grid
    using bin_container_type = dynamic_bin_container<bin_t, containers>;

    // Vecmem based view type
    using view_type =
        dmulti_view<dvector_view<bin_data_t>, dvector_view<entry_t>>;
    using const_view_type = dmulti_view<dvector_view<const bin_data_t>,
                                        dvector_view<const entry_t>>;

    // Vecmem based buffer type
#ifndef DETRAY_COMPILE_VITIS
    using buffer_type =
        dmulti_buffer<dvector_buffer<bin_data_t>, dvector_buffer<entry_t>>;
#endif

    /// Default constructor
    bin_storage() = default;
    /// Copy constructor
    bin_storage(const bin_storage&) = default;
    /// Move constructor
    bin_storage(bin_storage&&) = default;

    /// Construct containers using a memory resources
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST bin_storage(vecmem::memory_resource& resource)
        : m_bin_data(&resource), m_entry_data(&resource) {}
#endif // DETRAY_COMPILE_VITIS

    /// Construct grid data from containers - move
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type&& bin_data)
        : m_bin_data(std::move(bin_data.bins)),
          m_entry_data(std::move(bin_data.entries)) {}

    /// Construct the non-owning type from the @param offset into the global
    /// containers @param bin_data and the number of bins @param size
    template <bool owner = is_owning, std::enable_if_t<!owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type& bin_data, dindex offset,
                                   dindex size)
        : m_bin_data(bin_data.bins, dindex_range{offset, offset + size}),
          m_entry_data(
              bin_data.entries,
              dindex_range{0u, static_cast<dindex>(bin_data.entries.size())}) {}

    /// Construct bin storage from its vecmem view
    template <typename view_t,
              typename std::enable_if_t<detail::is_device_view<view_t>::value ,
                                        bool> = true>
    DETRAY_HOST_DEVICE bin_storage(const view_t& view)
        : m_bin_data(detray::detail::get<0>(view.m_view)),
          m_entry_data(detray::detail::get<1>(view.m_view)) {}

    /// Copy assignment
    bin_storage& operator=(const bin_storage&) = default;
    /// Move assignment
    bin_storage& operator=(bin_storage&&) = default;

    const bin_range_t& bin_data() const { return m_bin_data; }
    const entry_range_t& entry_data() const { return m_entry_data; }

    /// begin and end of the bin range
    /// @{
    DETRAY_HOST_DEVICE
    auto begin() {
        return iterator_adapter<bin_iterator_t>{
            detray::ranges::begin(m_bin_data), m_entry_data.data()};
    }
    DETRAY_HOST_DEVICE
    auto begin() const {
        return iterator_adapter<const_bin_iterator_t>{
            detray::ranges::cbegin(m_bin_data), m_entry_data.data()};
    }
    DETRAY_HOST_DEVICE
    auto end() {
        return iterator_adapter<bin_iterator_t>{detray::ranges::end(m_bin_data),
                                                m_entry_data.data()};
    }
    DETRAY_HOST_DEVICE
    auto end() const {
        return iterator_adapter<const_bin_iterator_t>{
            detray::ranges::cend(m_bin_data), m_entry_data.data()};
    }
    /// @}

    /// @returns the vecmem view of the bin storage
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return view_type{detray::get_data(m_bin_data),
                         detray::get_data(m_entry_data)};
    }
#endif // DETRAY_COMPILE_VITIS

    /// @returns the vecmem view of the bin storage - const
#ifndef DETRAY_COMPILE_VITIS
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return const_view_type{detray::get_data(m_bin_data),
                               detray::get_data(m_entry_data)};
    }
#endif // DETRAY_COMPILE_VITIS

    private:
    /// Container that holds all bin data when owning or a view into an
    /// externally owned container
    bin_range_t m_bin_data{};
    /// Container that holds all bin entries when owning or a view into an
    /// externally owned container
    entry_range_t m_entry_data{};
};

}  // namespace detray::detail
