/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// @brief base type for the storage element that wraps a bin content in the
///        grid backend storage.
///
/// The storage element also facilitates to provide a specific 'view' on the bin
/// content. This is usefull e.g. if the bin content is an iterable collection
/// of entries. The populators determine the behaviour of bin population, while
/// the entry type of an content is determined by the grid.
template <typename content_t, template <typename> class view_t,
          typename T = void>
struct bin_base {
    content_t m_content{};
    using content_type = content_t;
    using view_type = view_t<content_t>;
    using const_view_type = view_t<const content_t>;

    /// @returns the content that a bin stores
    DETRAY_HOST_DEVICE
    auto view() const -> const_view_type { return const_view_type(m_content); }

    /// Access to the bin content
    DETRAY_HOST_DEVICE
    auto view() -> view_type { return view_type(m_content); }

    /// @returns the bin content - const
    DETRAY_HOST_DEVICE
    auto content() const -> const content_t & { return m_content; }

    /// @returns the bin content
    DETRAY_HOST_DEVICE
    auto content() -> content_t & { return m_content; }

    /// Insert a @param new_content into the backend storage - move
    DETRAY_HOST_DEVICE
    void insert(content_t &&new_content) { m_content = std::move(new_content); }

    /// Insert a @param new_content into the backend storage - copy
    DETRAY_HOST_DEVICE
    void insert(const content_t &new_content) { m_content = new_content; }
};

/// @brief default implementation of an element in the grid backend storage
template <typename content_t, template <typename> class view_t>
struct bin
    : public bin_base<
          content_t, view_t,
          typename std::enable_if_t<std::is_standard_layout_v<content_t>>> {};

/// A replace populator that overrides whatever current content is in the bin
/// with a new entry.
///
/// @note entry type and bin content type are identicial in this case
struct replacer {

    /// Sort the entries contained in a bin content when fetched: not needed,
    /// as there is only a single entry
    static constexpr bool do_sort = false;

    template <typename content_t>
    using bin_type = bin<content_t, single_iterator>;

    /// Replace the bin content with a new entry - forwarding
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin index to be populated
    /// @param content new content to be added
    template <typename bin_storage_t, typename content_t>
    DETRAY_HOST_DEVICE void operator()(bin_storage_t &storage,
                                       const dindex gbin,
                                       content_t &&entry) const {
        storage[gbin].insert(std::forward<content_t>(entry));
    }

    /// Fetch a bin content from a storage element in the backend storage
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin index to be fetched
    ///
    /// @return a const iterator view on the bin content
    template <typename bin_storage_t>
    DETRAY_HOST_DEVICE auto view(const bin_storage_t &storage,
                                 const dindex gbin) const {
        return storage[gbin].view();
    }

    /// @returns an initialized bin in the backend storage
    template <typename content_t>
    DETRAY_HOST_DEVICE static constexpr auto init(
        content_t content = detail::invalid_value<content_t>())
        -> bin_type<content_t> {
        return {{content}};
    }
};

/// A complete populator that adds elements to a bin content which contains an
/// array of entries, until it is completed - ignored afterwards.
///
/// @tparam kDIM the dimension of the array in the bin content
/// @tparam entry_t the type of a single element in the bin content
/// @tparam kSORT sort the entries in the bin content
///
/// @note the entry type must be default constructible, since its default
/// entry may be interpreted as the empty element in which to put the new entry.
template <std::size_t kDIM = 1, bool kSORT = false,
          template <typename, std::size_t> class array_t = darray>
struct completer {

    /// Sort the entries contained in a bin content when viewed
    static constexpr bool do_sort = kSORT;

    template <typename entry_t>
    using bin_type = bin<array_t<entry_t, kDIM>, standard_iterator>;

    /// Complete the bin content with a new entry - copy
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin index to be populated
    /// @param entry new entry to complete the bin array with
    template <typename bin_storage_t, typename entry_t>
    DETRAY_HOST_DEVICE void operator()(bin_storage_t &storage,
                                       const dindex gbin,
                                       const entry_t entry) const {
        for (entry_t &stored_entry : storage[gbin].content()) {
            if (stored_entry == detail::invalid_value<entry_t>()) {
                stored_entry = entry;
            }
        }
    }

    /// Fetch a bin content from a storage element in the backend storage
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin to be viewed
    ///
    /// @return a const iterator view on the bin content
    template <typename bin_storage_t>
    DETRAY_HOST_DEVICE decltype(auto) view(const bin_storage_t &storage,
                                           const dindex gbin) const {
        return storage[gbin].view();
    }

    /// @returns an initialized storage element that contains @param entry in
    /// the first entry and sets all other entries to 'invalid'
    template <typename entry_t>
    DETRAY_HOST_DEVICE static constexpr auto init(
        entry_t entry = detail::invalid_value<entry_t>()) -> bin_type<entry_t> {
        // Initialize the storage element
        array_t<entry_t, kDIM> stored{};
        std::fill(stored.begin(), stored.end(),
                  detail::invalid_value<entry_t>());

        // The first entry always exists
        stored[0] = entry;

        return {{stored}};
    }
};

/// A regular attach populator that adds a new entry to a given bin.
///
/// @tparam kDIM the dimension of the bin content (array)
/// @tparam entry_t the type of a single element in the bin content
/// @tparam kSORT sort the entries in the bin
///
/// @note the number of entries are assumed to be the same in every bin.
/// @note the entry type must be default constructible, since its default
/// entry may be interpreted as the empty element in which to put the new entry.
template <std::size_t kDIM = 1, bool kSORT = false,
          template <typename, std::size_t> class array_t = darray>
struct regular_attacher {

    /// Sort the entries contained in a bin entry when viewed
    static constexpr bool do_sort = kSORT;

    template <typename entry_t>
    using bin_type = bin<array_t<entry_t, kDIM>, standard_iterator>;

    /// Append a new entry to the bin - forwarding
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin index to be populated
    /// @param entry new entry to add to the bin
    template <typename bin_storage_t, typename entry_t>
    DETRAY_HOST_DEVICE void operator()(bin_storage_t &storage,
                                       const dindex gbin,
                                       entry_t &&entry) const {
        for (auto &stored : storage[gbin].content()) {
            if (stored ==
                detail::invalid_value<std::remove_reference_t<entry_t>>()) {
                stored = std::forward<entry_t>(entry);
                break;
            }
        }
    }

    /// Fetch a bin from a storage element in the backend storage - const
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin index to be viewed
    ///
    /// @return a const iterator view on the bin content
    template <typename bin_storage_t>
    DETRAY_HOST_DEVICE auto view(const bin_storage_t &storage,
                                 const dindex gbin) const {
        return storage[gbin].view();
    }

    /// @returns an initialized storage element that contains @param entry in
    /// the first element of the bin
    template <typename entry_t>
    DETRAY_HOST_DEVICE static constexpr auto init(
        entry_t entry = detail::invalid_value<entry_t>()) -> bin_type<entry_t> {
        // Initialize the storage element
        array_t<entry_t, kDIM> stored{};
        std::fill(stored.begin(), stored.end(),
                  detail::invalid_value<entry_t>());

        // The first element always exists
        stored[0] = entry;

        return {{stored}};
    }
};

/// An irregular attach populator that adds the new entry to the collection of
/// entries in a given bin: each bin is dynamically sized.
///
/// @tparam kDIM the dimension of the array in the bin content
/// @tparam entry_t the type of a single element in the bin content
/// @tparam kSORT sort the elements in the bin content
///
/// @note since this type performs a bin search upon insertion and bin lookup,
/// it is comparatively slow and should only be used if there is a large
/// difference in the number of elements in the bin entries. Otherwise, use the
/// @c regular_attacher.
/// @note the entry type must be default constructible, since its default
/// entry may be interpreted as the empty element in which to put the new entry.
/*template <typename entry_t = dindex, std::size_t kDIM = 1, bool kSORT = false,
          template <typename...> class vector_t = dvector>
struct irregular_attacher {

    /// The storage element saves additional information to find the correct
    /// bin in a 'jagged' memory layout.
    template <typename T, template <typename> class V>
    struct indexed_bin
        : public bin_base<
              T, V, typename std::enable_if_t<std::is_standard_layout_v<T>>> {
        using base_type = bin_base<T, V, void>;
        using view_type = typename base_type::view_type;
        using const_view_type = typename base_type::view_type;
        /// The global bin a single entry in a bin content is associated with
        dindex gbin;

        /// Set to invalid entries by default
        indexed_bin() : gbin{dindex_invalid} {
            this->entry = detail::invalid_value<entry_t>();
        }
    };

    /// Sort the entries contained in a bin entry when fetched
    static constexpr bool do_sort = kSORT;

    template <typename entry_t>
    using bin_type = bin<vector_t<entry_t, kDIM>, single_iterator>;

    /// Add a new entry to the storage element - move
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin to be populated
    /// @param entry new entry to be added to the vector
    template<typename bin_storage_t, typename entry_t>
    DETRAY_HOST
    void operator()(bin_storage_t &storage,
                    const dindex &gbin, entry_t &&entry) const {
        auto last_gbin_entry = detail::upper_bound(storage.begin(),
storage.end(), entry); *last_gbin_entry = std::move(entry);
    }

    /// Add a new entry to the storage element - copy
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin to be populated
    /// @param entry new entry to be added to the vector
    template<typename bin_storage_t, typename entry_t>
    DETRAY_HOST
    void operator()(bin_storage_t &storage,
                    const dindex gbin, entry_t &entry) const {
        auto last_gbin_entry = detail::upper_bound(storage.begin(),
storage.end(), entry); *last_gbin_entry = entry;
    }

    /// Fetch a bin content from a storage element in the backend storage -
const
    ///
    /// @param storage the grid backend storage
    /// @param gbin the global grid bin to be fetched
    /// @param sub_idx fetch a single entry from the bin content
    ///
    /// @return the sub-content in the global bin
    template<typename bin_storage_t>
    DETRAY_HOST_DEVICE
    auto view(const bin_storage_t &storage,
               const dindex gbin) const { auto first_gbin_entry =
detail::lower_bound(storage.begin(), storage.end(), entry); auto last_gbin_entry
= detail::upper_bound(storage.begin(), storage.end(), entry);

    }

    /// @returns an initialized bin entry
    template<typename entry_t>
    DETRAY_HOST_DEVICE
    static constexpr auto init(
        entry_t content = detail::invalid_value<entry_t>()) -> bin_type<entry_t>
{
        // Set capacity the storage element
        entry_t stored(nentries);
        std::fill(stored.begin(), stored.end(),
detail::invalid_value<entry_t>());

        stored.at(0) = entry;

        return {stored};
    }
};*/

}  // namespace detray
