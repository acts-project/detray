/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s).
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

/// Does nothing
struct brute_force_view : public detail::dbase_view {
    brute_force_view() = default;
    brute_force_view(
        std::_Head_base<0, detray::brute_force_view, false> /*base*/) {}

    bool m_view = true;
};

/// @brief A surface finder that returns all surfaces in a volume (brute force)
struct brute_force_finder {

    using size_type = dindex;
    using view_type = brute_force_view;
    using const_view_type = brute_force_view;

    /// Default constructor
    constexpr brute_force_finder() = default;

    /// Constructor from memory resource: Not needed
    DETRAY_HOST
    constexpr brute_force_finder(vecmem::memory_resource * /*mr*/) {}

    /// Constructor from a vecmem view: Not needed
    DETRAY_HOST_DEVICE
    constexpr brute_force_finder(const brute_force_view & /*view*/) {}

    /// @returns the complete surface range of the search volume
    template <typename detector_t, typename track_t>
    DETRAY_HOST_DEVICE constexpr auto search(
        const detector_t & /*det*/,
        const typename detector_t::volume_type &volume,
        const track_t & /*track*/) const -> dindex_range {
        return volume.full_range();
    }

    /// @return the view on the brute force finder
    constexpr auto get_data() const noexcept -> brute_force_view { return {}; }

    /// @note outside of navigation, the number of elements is unknown
    constexpr auto size() const noexcept -> std::size_t { return 0; }

    /// @note outside of navigation, the number of elements is unknown
    constexpr auto empty() const noexcept -> bool { return true; }

    /// Does nothing
    template <typename... Args>
    constexpr auto push_back(Args &&... /*args*/) const noexcept -> void {}
};

}  // namespace detray