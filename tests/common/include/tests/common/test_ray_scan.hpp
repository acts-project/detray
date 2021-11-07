/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "detray/utils/ray_gun.hpp"

namespace detray {

/** Check if a set of volume index pairs form a trace.
 *
 * @tparam the record entry, which must contain the portal index and link
 *
 * @param volume_records the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @note Empty traces pass automatically.
 *
 * @return true if the volumes form a connected chain.
 */
template <typename entry_type = std::pair<dindex, dindex>>
inline bool check_connectivity(
    std::vector<std::pair<entry_type, entry_type>> volume_records,
    dindex start_volume = 0) {
    // Keep record of leftovers
    std::stringstream record_stream;

    // Where are we on the trace?
    dindex on_volume = start_volume;

    // Init chain search
    // find connection record for the current volume
    auto record =
        find_if(volume_records.begin(), volume_records.end(),
                [&](const std::pair<entry_type, entry_type> &rec) -> bool {
                    return (std::get<1>(rec.first) == on_volume) or
                           (std::get<1>(rec.second) == on_volume);
                });

    while (record != volume_records.end()) {
        auto first_vol = std::get<1>(record->first);
        auto second_vol = std::get<1>(record->second);

        record_stream << "On volume: " << on_volume << " and record ("
                      << first_vol << ", " << second_vol << ")";

        // update to next volume
        on_volume = on_volume == first_vol ? second_vol : first_vol;

        record_stream << "-> next volume: " << on_volume << std::endl;

        // Don't search this key again -> only one potential key with current
        // index left
        volume_records.erase(record);

        // find connection record for the current volume
        record =
            find_if(volume_records.begin(), volume_records.end(),
                    [&](const std::pair<entry_type, entry_type> &rec) -> bool {
                        return (std::get<1>(rec.first) == on_volume) or
                               (std::get<1>(rec.second) == on_volume);
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
    // obj id and obj mother volume
    using trace_entry = std::pair<dindex, dindex>;
    // Indices of volumes that are linked by portals
    std::vector<std::pair<trace_entry, trace_entry>> portal_trace = {};
    // Indices of volumes per module surface
    std::vector<trace_entry> surface_trace = {};
    // Debug output if an error in the record is discovered
    std::stringstream record_stream;

    // Always read 2 elements from the sorted records vector
    struct record_doublet {
        // Two times {portal index, portal intersection}
        const typename record_type::value_type &lower, upper;

        // Is the record doublet we picked up made up if a portal and a surface?
        inline bool is_mixed() const {
            return (lower.second.index != lower.second.link and
                    upper.second.index == upper.second.link) or
                   (lower.second.index == lower.second.link and
                    upper.second.index != upper.second.link);
        }

        // Is this doublet connected via a portal intersection? (the second
        // requirement guarantees that this indeed a portal crossing, i.e.
        // changeing volumes)
        inline bool operator()() const {
            return (lower.second == upper.second) and
                   (lower.second.index != upper.second.index);
        }
    };

    auto is_portal = [](const std::pair<dindex, intersection> &record) -> bool {
        // A portal links to another volume thab it belongs to
        return record.second.index != record.second.link;
    };

    // Don't look at the end volume, as it is not connected at one side
    for (size_t rec = 0; rec < volume_record.size() - 1;) {

        if (not is_portal(volume_record.at(rec))) {
            surface_trace.emplace_back(volume_record[rec].first,
                                       volume_record[rec].second.index);
            rec++;
            continue;
        }

        // Get 2 possibly connected entries
        record_doublet doublet = {.lower = volume_record.at(rec),
                                  .upper = volume_record.at(rec + 1)};

        // Advance to inspect next pair
        rec += 2;

        record_stream << doublet.lower.second.index << " ("
                      << doublet.lower.second.path << ")" << std::endl;
        record_stream << doublet.upper.second.index << " ("
                      << doublet.upper.second.path << ")" << std::endl;

        if (doublet()) {
            // Insert into set of edges
            trace_entry lower{doublet.lower.first, doublet.lower.second.index};
            trace_entry upper{doublet.upper.first, doublet.upper.second.index};
            portal_trace.emplace_back(lower, upper);
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

            return std::make_pair(portal_trace, surface_trace);
        }
    }

    // Look at the last entry
    if (not is_portal(volume_record.back())) {
        std::cerr << "We don't leave the detector by portal!" << std::endl;
    }
    // Put this in the surface trace, because it is not part of a doublet
    else {
        trace_entry lower{volume_record.back().first,
                          volume_record.back().second.index};
        trace_entry upper{volume_record.back().first,
                          volume_record.back().second.link};
        portal_trace.emplace_back(lower, upper);
    }

    return std::make_pair(portal_trace, surface_trace);
}

/** Build an adjacency list from volume traces.
 *
 * @tparam record_type container that contains volume idx, intersection pairs
 *
 * @param volume_record the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @return a set of volume connections that were found by portal intersection
 *         of a ray.
 */
template <typename portal_record, typename surface_record,
          typename entry_type = std::pair<dindex, dindex>,
          std::enable_if_t<std::is_same_v<typename portal_record::value_type,
                                          std::pair<entry_type, entry_type>>,
                           bool> = true,
          std::enable_if_t<
              std::is_same_v<typename surface_record::value_type, entry_type>,
              bool> = true>
inline auto build_adjacency(
    const portal_record &portal_trace, const surface_record &surface_trace,
    std::map<dindex, std::map<dindex, dindex>> &adj_list = {},
    std::unordered_set<dindex> &obj_hashes = {}) {

    // Every surface that was recorded adds a link to the mother volume
    for (const auto &record : surface_trace) {
        const auto sf_index = std::get<0>(record);
        const auto vol_index = std::get<1>(record);
        // Check whether we have seen this surface in this volume before
        if (obj_hashes.find(sf_index) == obj_hashes.end()) {
            adj_list[vol_index][vol_index]++;
            obj_hashes.insert(sf_index);
        }
    }

    // Portal in first volume links to second volume in the record
    for (const auto &record : portal_trace) {
        const auto pt_index_1 = std::get<0>(record.first);
        const auto vol_index_1 = std::get<1>(record.first);
        const auto pt_index_2 = std::get<0>(record.second);
        const auto vol_index_2 = std::get<1>(record.second);

        if (obj_hashes.find(pt_index_1) == obj_hashes.end()) {
            adj_list[vol_index_1][vol_index_2]++;
            obj_hashes.insert(pt_index_1);
        }
        // Assume the return link for now (filter out portal that leaves world)
        if (vol_index_2 != dindex_invalid) {
            if (obj_hashes.find(pt_index_2) == obj_hashes.end()) {
                adj_list[vol_index_2][vol_index_1]++;
                obj_hashes.insert(pt_index_2);
            }
        }
    }

    return adj_list;
}

}  // namespace detray