/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <algorithm>
#include <cmath>
#include <iterator>
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

/** Check if a set of volume indices from portal intersections from a path
 *  (works even if the pairs are not sorted).
 *
 * @tparam the record entry, which must contain the portal and mother volume
 *         index
 *
 * @param volume_records the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @return true if the volumes indices form a connected chain.
 */
template <bool check_sorted_trace = true,
          typename entry_type = std::pair<dindex, dindex>>
inline bool check_connectivity(
    std::vector<std::pair<entry_type, entry_type>> volume_records,
    dindex start_volume = 0) {
    // Keep record of leftovers
    std::stringstream record_stream;

    // Where are we on the trace?
    dindex on_volume = start_volume;

    // If the intersection trace comes from the ray gun/ trace intersections
    // function it should be sorted, which is the stronger constraint
    using records_iterator_t = decltype(volume_records.begin());
    std::function<records_iterator_t(dindex)> get_record;
    if constexpr (check_sorted_trace) {
        // Get the next record
        get_record = [&](dindex next) -> records_iterator_t {
            return volume_records.begin() + next;
        };
    } else {
        // Search for the existence of a fitting record
        get_record = [&](dindex /*next*/) -> records_iterator_t {
            return find_if(
                volume_records.begin(), volume_records.end(),
                [&](const std::pair<entry_type, entry_type> &rec) -> bool {
                    return (std::get<1>(rec.first) == on_volume) or
                           (std::get<1>(rec.second) == on_volume);
                });
        };
    }

    // Init chain search
    dindex i = 0;
    auto record = get_record(i);

    // Check first volume index, which has no partner otherwise
    if ((not std::get<1>(record->first) == start_volume)) {
        std::cerr << "First record does not start at given initial volume: "
                  << std::get<1>(record->first) << " vs. " << start_volume
                  << std::endl;

        return false;
    }

    // Walk along the trace
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
        if constexpr (not check_sorted_trace) {
            volume_records.erase(record);
        }

        // find connection record for the current volume
        record = get_record(++i);
    }

    // There are unconnected elements left
    if (not on_volume == dindex_invalid) {
        std::cerr << "In trace finding: " << std::endl;
        std::cerr << record_stream.str() << std::endl;

        return false;
    }

    return true;
}

/** Check if a recording of portal/surface intersections are form a coherent
 *  trace through a geometry. The linking information the intersection holds
 *  is only used to sort surfaces from portals.
 *
 * @tparam record_type container that contains object indices and intersections
 *
 * @param volume_record the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @note the input record needs to be sorted according to the distance from the
 *       ray origin
 *
 * @return a set of volume connections that were found by portal intersection
 *         of a ray.
 */
template <typename record_container = dvector<std::pair<dindex, intersection>>>
inline auto trace_intersections(const record_container &volume_record,
                                dindex start_volume = 0) {
    // obj id and obj mother volume
    using trace_entry = std::pair<dindex, dindex>;
    // Indices of volumes that are linked by portals
    std::vector<std::pair<trace_entry, trace_entry>> portal_trace = {};
    // Indices of volumes per module surface
    std::vector<trace_entry> surface_trace = {};
    // Debug output if an error in the record is discovered
    std::stringstream record_stream;

    // Readable access to the data of a recorded intersection
    struct record {
        using type = typename record_container::value_type;
        const type &entry;

        // getter
        inline auto &object_id() const { return entry.first; }
        inline auto &inters() const { return entry.second; }
        inline auto &volume_id() const { return entry.second.index; }
        inline auto &volume_link() const { return entry.second.link; }
        inline auto &dist() const { return entry.second.path; }

        inline bool is_portal() const {
            // A portal links to another volume than it belongs to
            return entry.second.index != entry.second.link;
        }

        inline bool is_portal(
            const std::pair<dindex, intersection> &inters_pair) const {
            // A portal links to another volume than it belongs to
            const record rec{inters_pair};
            return rec.volume_id() != rec.volume_link();
        }
    };

    for (size_t rec = 0; rec < volume_record.size() - 1;) {

        // For portals read 2 elements from the sorted records vector
        const record current_rec = record{volume_record.at(rec)};
        const record next_rec = record{volume_record.at(rec + 1)};

        // Add an entry to the surface trace and continue more fine-grained
        if (not current_rec.is_portal()) {
            surface_trace.emplace_back(current_rec.object_id(),
                                       current_rec.volume_id());
            rec++;
            continue;
        }

        record_stream << current_rec.volume_id() << " (" << current_rec.dist()
                      << ")" << std::endl;
        record_stream << next_rec.volume_id() << " (" << next_rec.dist() << ")"
                      << std::endl;

        // Is this doublet connected via a portal intersection? (the second
        // requirement guarantees that this indeed a portal crossing, i.e.
        // changeing volumes)
        const bool is_valid = (current_rec.inters() == next_rec.inters()) and
                              (current_rec.volume_id() != next_rec.volume_id());

        // Is the record doublet we picked up made up of a portal and a surface?
        const bool is_mixed =
            (current_rec.is_portal() and not next_rec.is_portal()) or
            (next_rec.is_portal() and not current_rec.is_portal());

        if (is_valid and not is_mixed) {
            // Insert into set of edges
            trace_entry lower{current_rec.object_id(), current_rec.volume_id()};
            trace_entry upper{next_rec.object_id(), next_rec.volume_id()};
            portal_trace.emplace_back(lower, upper);
        }
        // Something went wrong
        else {
            // Print search log
            std::cerr << "<<<<<<<<<<<<<<<" << std::endl;
            std::cerr << "volume idx (distance from ray origin)\n" << std::endl;

            std::cerr << record_stream.str() << std::endl;

            std::cerr << "Ray terminated at portal x-ing " << (rec + 1) / 2
                      << ": (" << current_rec.object_id() << ", "
                      << current_rec.dist() << ") <-> (" << next_rec.object_id()
                      << ", " << next_rec.dist() << ")" << std::endl;
            std::cerr << std::boolalpha << "Portals are connected: " << is_valid
                      << std::endl;
            std::cerr << std::boolalpha
                      << "Portals link to themselves: " << is_mixed
                      << std::endl;

            record rec_front{volume_record.front()};
            record rec_back{volume_record.back()};
            std::cerr << "Start volume : " << start_volume << std::endl;
            std::cerr << "- first intersection: (" << rec_front.object_id()
                      << ", " << rec_front.dist() << ")," << std::endl;
            std::cerr << "- last intersection: (" << rec_back.object_id()
                      << ", " << rec_back.dist() << ")," << std::endl;

            std::cerr << ">>>>>>>>>>>>>>>" << std::endl;

            return std::make_pair(portal_trace, surface_trace);
        }

        // Advance to inspect next pair
        rec += 2;
    }

    record rec_back{volume_record.back()};
    // Look at the last entry
    if (not rec_back.is_portal()) {
        std::cerr << "We don't leave the detector by portal!" << std::endl;
    } else {
        trace_entry lower(rec_back.object_id(), rec_back.volume_id());
        trace_entry upper(rec_back.object_id(), rec_back.volume_link());
        portal_trace.emplace_back(lower, upper);
    }

    return std::make_pair(portal_trace, surface_trace);
}

/** Build an adjacency list from intersection traces.
 *
 * @tparam portal_trace_type container of portal link pairs
 * @tparam surface_trace_type container of surface links
 *
 * @param portal_trace the portal indices and their volume links (in adjacent
 *                     portal pairs)
 * @param surface_trace the surface indices and their volume links
 * @param obj_hashes record which surfaces/portals were already added
 *
 * @return an adjacency list from the traced ray scan of a given geometry.
 */
template <
    typename portal_trace_type, typename surface_trace_type,
    typename entry_type = std::pair<dindex, dindex>,
    std::enable_if_t<std::is_same_v<typename portal_trace_type::value_type,
                                    std::pair<entry_type, entry_type>>,
                     bool> = true,
    std::enable_if_t<
        std::is_same_v<typename surface_trace_type::value_type, entry_type>,
        bool> = true>
inline auto build_adjacency(
    const portal_trace_type &portal_trace,
    const surface_trace_type &surface_trace,
    std::map<dindex, std::map<dindex, dindex>> &adj_list,
    std::unordered_set<dindex> &obj_hashes) {

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