/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

/// @brief A surface finder that returns all surfaces in a volume (brute force)
struct brute_force_finder {

    /// Default constructor
    brute_force_finder() = default;

    /// Constructor from memory resource: Not needed
    DETRAY_HOST
    brute_force_finder(vecmem::memory_resource & /*mr*/) {}

    /// Constructor from a vecmem view: Not needed
    template <typename view_t,
              std::enable_if_t<
                  !std::is_same_v<brute_force_finder, view_t> &&
                      !std::is_base_of_v<vecmem::memory_resource, view_t>,
                  bool> = true>
    DETRAY_HOST_DEVICE brute_force_finder(const view_t & /*view*/) {}

    /// @returns the complete surface range of the search volume
    template <typename detector_t, typename track_t>
    DETRAY_HOST_DEVICE dindex_range
    search(const detector_t & /*det*/,
           const typename detector_t::volume_type &volume,
           const track_t & /*track*/) const {
        return volume.range();
    }
};

}  // namespace detray
