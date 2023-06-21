/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algorithms.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/detail/populator_impl.hpp"

namespace detray {

/// Enforce the populator interface and implement some common functionality
/// @todo remove this interface
template <typename populator_impl_t>
class populator {

    public:
    // Insert the correct container/algebra types
    using impl = populator_impl_t;

    template <typename entry_t>
    using bin_type = typename impl::template bin_type<entry_t>;

    /// Populate bin with a new entry - forwarding
    ///
    /// @param storage the global bin storage
    /// @param gbin the global bin index to be populated
    /// @param entry the single entry to be added to the bin content
    template <typename serialized_storage, typename entry_t>
    DETRAY_HOST_DEVICE void operator()(serialized_storage &storage,
                                       const dindex gbin,
                                       entry_t &&entry) const {
        m_populator_impl(storage, gbin, std::forward<entry_t>(entry));

        if constexpr (impl::do_sort) {
            detail::sequential_sort(storage[gbin].content().begin(),
                                    storage[gbin].content().end());
        }
    }

    /// Fetch a bin entry from the grid backend storage
    template <typename serialized_storage>
    DETRAY_HOST_DEVICE constexpr auto view(const serialized_storage &storage,
                                           const dindex gbin) const {
        return m_populator_impl.view(storage, gbin);
    }

    /// Fetch a bin entry from the grid backend storage
    template <typename serialized_storage>
    DETRAY_HOST_DEVICE auto view(serialized_storage &storage,
                                 const dindex gbin) {
        return m_populator_impl.view(storage, gbin);
    }

    /// @returns a default initialized bin entry
    template <typename entry_t>
    DETRAY_HOST_DEVICE static constexpr auto init() {
        return impl::template init<entry_t>();
    }

    /// @returns an initialized bin entry with @param entry.
    template <typename entry_t>
    DETRAY_HOST_DEVICE static constexpr auto init(entry_t entry) {
        return impl::init(entry);
    }

    impl m_populator_impl{};
};

}  // namespace detray
