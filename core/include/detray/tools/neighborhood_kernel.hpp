/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <utility>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

struct neighborhood_kernel {

    /// Variadic unrolled neighborhood lookups - any integer sequence
    ///
    /// @tparam track_t The type of the track
    /// @tparam finder_container_t The type of the type of the finder container
    /// @tparam finder_range_t The surface finder range type
    /// @tparam first_finder_id The first surface finder group id
    ///
    /// @param tack the track information
    /// @param surface_finders the masks container
    /// @param finder_id the range within the mask group to be checked
    /// @param finder_range the current mask group id
    /// @param available_ids the finder ids to be checked (only needed to set
    ///                      the first id for the call)
    ///
    /// @return a collection of neighboring surface indices (empty if no
    /// surface finder could be matched)
    template <typename track_t, typename finder_container_t,
              typename finder_range_t, unsigned int first_finder_id,
              unsigned int... remaining_finder_ids>
    DETRAY_HOST_DEVICE inline auto unroll_search(
        const track_t &track, const finder_container_t &surface_finders,
        const typename finder_container_t::id_type finder_id,
        const finder_range_t &finder_range,
        std::integer_sequence<unsigned int, first_finder_id,
                              remaining_finder_ids...>
        /*available_ids*/) {

        // Pick the first one for interseciton
        if (finder_id == first_finder_id) {
            return;
        }

        // The reduced integer sequence
        std::integer_sequence<unsigned int, remaining_finder_ids...> remaining;

        // Unroll as long as you have at least 1 entries
        if constexpr (remaining.size() >= 1) {
            return (unroll_intersect(track, surface_finders, finder_id,
                                     finder_range, remaining));
        }

        // No intersection was found
        return;
    }

    /// Kernel method that does a surface neighborhood lookup given a surface
    ///  finder link provided by a volume.
    ///
    /// @tparam volume_t The calling volume
    /// @tparam track_t The type of the track
    /// @tparam sf_finder_container_t The type of the surface finder container
    ///
    /// @param track the track information
    /// @param volume the mother volume
    ///
    /// @return an collection of surface indices (invalid if no finder matches)
    template <typename volume_t, typename track_t, typename finder_container_t>
    inline auto operator()(const volume_t &volume,
                                 const track_t & /*track*/,
                                 const finder_container_t & /*sf_finders*/) {

        // Unroll the searching depending on the number of surface finder types
        // using finder_defs = typename volume_t::sf_finders;
        return volume.range();
        /*return unroll_search(
            track, sf_finders, finder_id, finder_range,
            std::make_integer_sequence<unsigned int, sf_finders::n_types>{});*/
    }
};
}  // namespace detray