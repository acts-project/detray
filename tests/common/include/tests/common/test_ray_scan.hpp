/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "utils/ray_gun.hpp"

namespace detray {

/** Check if a set of volume index pairs form a trace.
 *
 * @param volume_records the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @note Empty traces pass automatically.
 *
 * @return true if the volumes form a connected chain.
 */
inline bool check_connectivity(
    std::set<std::pair<dindex, dindex>> volume_records,
    dindex start_volume = 0) {
    // Keep record of leftovers
    std::stringstream record_stream;

    // Where are we on the trace?
    dindex on_volume = start_volume;

    // Init chain search
    // find connection record for the current volume
    auto record = find_if(volume_records.begin(), volume_records.end(),
                          [&](const std::pair<dindex, dindex> &rec) -> bool {
                              return (rec.first == on_volume) or
                                     (rec.second == on_volume);
                          });

    while (record != volume_records.end()) {

        record_stream << "On volume: " << on_volume << " and record ("
                      << record->first << ", " << record->second << ")";

        // update to next volume
        on_volume = on_volume == record->first ? record->second : record->first;

        record_stream << "-> next volume: " << on_volume << std::endl;

        // Don't search this key again -> only one potential key with current
        // index left
        volume_records.erase(record);

        // find connection record for the current volume
        record = find_if(volume_records.begin(), volume_records.end(),
                         [&](const std::pair<dindex, dindex> &rec) -> bool {
                             return (rec.first == on_volume) or
                                    (rec.second == on_volume);
                         });
    }

    // There are unconnected elements left
    if (not volume_records.empty()) {
        std::cerr << "In trace finding: " << std::endl;
        std::cerr << record_stream.str() << std::endl;

        return false;
    }

    return true;
}

/** Check if a vector of volume portal intersections are connected.
 *
 * @tparam record_type container that contains volume idx, intersection pairs
 *
 * @param volume_record the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @return a set of volume connections that were found by portal intersection
 *         of a ray.
 */
template <typename record_type = dvector<std::pair<dindex, intersection>>>
inline auto trace_volumes(const record_type &volume_record,
                          dindex start_volume = 0) {
    std::set<std::pair<dindex, dindex>> valid_volumes = {};
    std::stringstream record_stream;

    // Always read 2 elements from the sorted records vector
    struct record_doublet {
        // Two volume index, portal intersections
        const typename record_type::value_type &lower, upper;

        // Is this doublet connected via the portal intersection?
        inline bool operator()() const { return lower.second == upper.second; }
    };

    // Filter out intersection with non-portals
    for (auto rec : volume_record) {
    }

    // Excluding the volume where we leave world, the rest should come in pairs
    if ((volume_record.size() - 1) % 2) {
        return valid_volumes;
    }

    // Don't look at the end volume, as it is not connected at one side
    for (size_t rec = 0; rec < volume_record.size() - 1; rec += 2) {

        // Get 2 possibly connected entries
        record_doublet doublet = {.lower = volume_record[rec],
                                  .upper = volume_record[rec + 1]};

        record_stream << doublet.lower.first << " ("
                      << doublet.lower.second.path << ")" << std::endl;
        record_stream << doublet.upper.first << " ("
                      << doublet.upper.second.path << ")" << std::endl;

        if (doublet()) {
            // Insert into set of edges
            valid_volumes.emplace(doublet.lower.first, doublet.upper.first);
        }
        // Something went wrong
        else {
            // Print search log
            std::cerr << "<<<<<<<<<<<<<<<" << std::endl;
            std::cerr << "volume idx (distance from ray origin)\n" << std::endl;

            std::cerr << record_stream.str() << std::endl;

            std::cerr << "Ray terminated at portal x-ing " << (rec + 1) / 2
                      << ": (" << doublet.lower.first << ", "
                      << doublet.lower.second.path << ") <-> ("
                      << doublet.upper.first << ", "
                      << doublet.upper.second.path << ")" << std::endl;

            std::cerr << "Start volume : " << start_volume << std::endl;
            std::cerr << "- first intersection: ("
                      << volume_record.front().first << ", "
                      << volume_record.front().second.path << ")," << std::endl;
            std::cerr << "- last intersection: (" << volume_record.back().first
                      << ", " << volume_record.back().second.path << "),"
                      << std::endl;

            std::cerr << ">>>>>>>>>>>>>>>" << std::endl;

            return valid_volumes;
        }
    }

    return valid_volumes;
}

}  // namespace detray