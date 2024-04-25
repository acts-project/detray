/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace detray::detector_scanner {

/// Check if a set of volume indices from portal intersections from a path
/// (works even if the pairs are not sorted). That the volume links of the
/// portals at a single boundary crossing match was already checked by the
/// @c trace_intersections function at this point.
///
/// For example, with (a, b) modeling a valid portal crossing between volume 'a'
/// and 'b' (leaving out the portal indices that are also in those pairs), then:
///     (0, 1)(1, 16)(16, 17)  is valid, (0, 1)(1, 15)(16, 17) is not
///         |__|   |__|    |_                |__|   | x |   |_
///
/// Also checks that the first portal that was found, lies in the start volume.
///
/// @tparam invalid_value to keep the implementation simple, all indices are of
///                       type @c dindex here, but that does not need to be the
///                       case in the intersection code, so we need to know what
///                       the invalid value is interpreted as a @c dindex
/// @tparam check_sorted_trace if the trace is not sorted, perform a search for
///                            a fitting record accross the trace.
/// @tparam entry_type the record entry, which must contain the portal index
///                    and the volume index the portal was discovred in.
///
/// @param trace the recorded portal crossings between volumes
/// @param start_volume where the ray started
///
/// @return true if the volumes indices form a connected chain.
template <dindex invalid_value = dindex_invalid, bool check_sorted_trace = true,
          typename entry_type = std::pair<dindex, dindex>>
inline bool check_connectivity(
    std::vector<std::pair<entry_type, entry_type>> trace,
    dindex start_volume = 0u) {

    /// Error messages
    std::stringstream err_stream;

    /// Print errors of this function
    auto print_err = [](const std::stringstream &stream) {
        std::cerr << "\n<<<<<<<<<<<<<<< ERROR in connectivity check\n"
                  << std::endl;
        std::cerr << stream.str() << std::endl;
        std::cerr << "\n>>>>>>>>>>>>>>>\n" << std::endl;
    };

    // There must always be portals!
    if (trace.empty()) {
        err_stream << "Trace empty!";
        print_err(err_stream);

        return false;
    }
    // Keep record of leftovers
    std::stringstream record_stream;

    // Where are we on the trace?
    dindex current_volume = start_volume;

    // If the intersection trace comes from the ray gun/trace intersections
    // function it should be sorted, which is the stronger constraint
    using records_iterator_t = decltype(trace.begin());
    using index_t = typename records_iterator_t::difference_type;
    std::function<records_iterator_t(index_t)> get_connected_record;
    if constexpr (check_sorted_trace) {
        // Get the next record
        get_connected_record = [&](index_t next) -> records_iterator_t {
            auto rec = trace.begin() + next;
            // Make sure that the record contains the volume that is currently
            // being checked for connectivity
            if (rec != trace.end() &&
                ((std::get<1>(rec->first) == current_volume) ||
                 (std::get<1>(rec->second) == current_volume))) {
                return rec;
            }
            return trace.end();
        };
    } else {
        // Search for the existence of a fitting record over the entire trace
        get_connected_record = [&](index_t /*next*/) -> records_iterator_t {
            return find_if(
                trace.begin(), trace.end(),
                [&](const std::pair<entry_type, entry_type> &rec) -> bool {
                    return (std::get<1>(rec.first) == current_volume) or
                           (std::get<1>(rec.second) == current_volume);
                });
        };
    }

    // Init chain search
    index_t i{0};
    auto record = get_connected_record(i);

    // Check first volume index, which has no partner otherwise
    if (std::get<1>(record->first) != start_volume) {
        err_stream << "First record does not start at given initial volume: "
                   << std::get<1>(record->first) << " vs. " << start_volume;

        print_err(err_stream);

        return false;
    }

    // Walk along the trace as long as a connection is found
    while (record != trace.end()) {
        auto first_vol = std::get<1>(record->first);
        auto second_vol = std::get<1>(record->second);

        record_stream << "On volume: " << current_volume << " and record ("
                      << first_vol << ", " << second_vol << ")";

        // update to next volume
        current_volume = (current_volume == first_vol ? second_vol : first_vol);

        record_stream << " -> next volume: " << current_volume << std::endl;

        // Don't search this key again -> only one potential key with current
        // index left
        if constexpr (not check_sorted_trace) {
            trace.erase(record);
        }

        // find connected record for the current volume
        record = get_connected_record(++i);
    }

    // There are unconnected elements left (we didn't leave world before
    // termination)
    if (current_volume != invalid_value) {
        err_stream
            << "Didn't leave world or unconnected elements left in trace:"
            << "\n\nValid connections that were found:" << std::endl;
        err_stream << record_stream.str();

        print_err(err_stream);

        return false;
    }

    return true;
}

/// Check if a recording of portal/module intersections form a coherent
/// trace through a geometry.
///
/// Various consistency checks are performed on the trace:
///     - Was a surface intersected as part of the wrong volume?
///     - Do the volume links at adjacent portals point at each other?
///     - Is a portal crossing happening within a volume, as opposed to at its
///       boundary surfaces?
///     - Do portals always appear in pairs (exit + entry portal)?
///
/// @tparam record_container contains volume indices and intersections
///
/// @param volume_record the recorded portal crossings between volumes
/// @param start_volume where the ray started
///
/// @note the input record needs to be sorted according to the distance from the
///       ray origin
///
/// @return a set of volume connections that were found by portal intersection
///         of a ray.
template <dindex invalid_value = dindex_invalid, typename record_container>
inline auto trace_intersections(const record_container &intersection_records,
                                dindex start_volume = 0u) {
    /// surface index and index of the volume the intersection was found in
    using trace_entry = std::pair<dindex, dindex>;
    /// Pairs of adjacent portals along the ray
    std::vector<std::pair<trace_entry, trace_entry>> portal_trace = {};
    /// Trace of module surfaces
    std::vector<trace_entry> module_trace = {};
    /// Debug output if an error in the trace is discovered
    std::stringstream record_stream;
    /// Error messages
    std::stringstream err_stream;
    bool error_code{true};

    /// Readable access to the data of a recorded intersection
    struct record {
        const typename record_container::value_type &entry;

        /// getter
        /// @{
        inline auto surface_idx() const {
            return entry.intersection.sf_desc.index();
        }
        inline auto surface_volume_idx() const {
            return entry.intersection.sf_desc.volume();
        }
        inline auto &inters() const { return entry.intersection; }
        inline auto &volume_idx() const { return entry.vol_idx; }
        inline auto &volume_link() const {
            return entry.intersection.volume_link;
        }
        inline auto &dist() const { return entry.intersection.path; }
        inline bool is_portal() const {
            return entry.intersection.sf_desc.is_portal();
        }
        inline bool is_sensitive() const {
            return entry.intersection.sf_desc.is_sensitive();
        }
        inline bool is_passive() const {
            return entry.intersection.sf_desc.is_passive();
        }
        /// @}
    };

    /// Print errors of this function
    auto print_err = [&error_code](const std::stringstream &stream) {
        std::cerr << "\n<<<<<<<<<<<<<<< ERROR intersection trace\n"
                  << std::endl;
        std::cerr << stream.str() << std::endl;
        std::cerr << "\n>>>>>>>>>>>>>>>\n" << std::endl;
        error_code = false;
    };

    // No intersections found by ray
    if (intersection_records.empty()) {
        err_stream << "No surfaces found!";
        print_err(err_stream);

        return std::make_tuple(portal_trace, module_trace, error_code);
    }

    // If there is only one surface in the trace, it must be a portal
    if (intersection_records.size() == 1u) {

        const record rec{intersection_records.at(0u)};

        // No exit potal
        if (not rec.is_portal()) {
            const std::string sf_type{rec.is_sensitive() ? "sensitive"
                                                         : "passive"};

            err_stream << "We don't leave the detector by portal!" << std::endl;
            err_stream << "Only found single " << sf_type
                       << " surface: portal(s) missing!";

            print_err(err_stream);

            return std::make_tuple(portal_trace, module_trace, error_code);
        }
    }

    // Go through recorded intersection (two at a time)
    dindex current_vol = start_volume;
    for (std::size_t rec = 0u; rec < (intersection_records.size() - 1u);) {

        const record current_rec = record{intersection_records.at(rec)};
        const record next_rec = record{intersection_records.at(rec + 1u)};

        // Keep a debug stream
        record_stream << current_rec.volume_idx() << "\t"
                      << current_rec.inters() << std::endl;

        // If the current record is not a portal, add an entry to the module
        // trace and continue in more fine-grained steps (sensitive/passive
        // surfaces do not come in pairs)
        if (not current_rec.is_portal()) {
            // Check that the surface was found in the volume it claims to
            // belong to
            const bool is_in_volume =
                (current_rec.volume_idx() ==
                 current_rec.surface_volume_idx()) and
                (current_rec.surface_volume_idx() == current_vol);
            if (is_in_volume) {
                module_trace.emplace_back(current_rec.surface_idx(),
                                          current_rec.volume_idx());
            } else {
                err_stream << "\n(!!) Surface " << current_rec.surface_idx()
                           << " outside of its volume (Found in: "
                           << current_vol << ", belongs in: "
                           << current_rec.surface_volume_idx() << ")\n";

                err_stream << record_stream.str();

                print_err(err_stream);

                return std::make_tuple(portal_trace, module_trace, error_code);
            }
            ++rec;
            continue;
        }
        // If the record is a portal, the current volume switches
        // the portals in the pair may be unordered, so check both
        else if (current_vol == current_rec.volume_idx()) {
            current_vol = current_rec.volume_link();
        } else if (current_vol == next_rec.volume_idx()) {
            current_vol = next_rec.volume_link();
        }

        record_stream << next_rec.volume_idx() << "\t" << next_rec.inters()
                      << std::endl;

        // Check that also the second surface was found in the volume it claims
        // to belong to
        const bool is_in_volume =
            next_rec.volume_idx() == next_rec.surface_volume_idx();
        // Is this doublet connected via a valid portal intersection?
        const bool is_valid =
            (current_rec.inters() == next_rec.inters()) and
            (current_rec.volume_idx() == next_rec.volume_link()) and
            (next_rec.volume_idx() == current_rec.volume_link());
        // Is this indeed a portal crossing, i.e. changing volumes?
        const bool is_self_link =
            current_rec.volume_idx() == next_rec.volume_idx();
        // Is the record doublet we picked made up of a portal and a module?
        const bool is_mixed =
            (current_rec.is_portal() and not next_rec.is_portal()) or
            (next_rec.is_portal() and not current_rec.is_portal());

        if (not is_in_volume) {
            record_stream << "\n(!!) Surface outside of its volume (Found: "
                          << next_rec.volume_idx()
                          << ", belongs in: " << next_rec.surface_volume_idx()
                          << ")" << std::endl;
        }
        if (not is_valid) {
            record_stream << "\n(!!) Not a valid portal crossing ("
                          << current_rec.volume_idx() << " <-> "
                          << next_rec.volume_idx() << "):\nPortals are not "
                          << "connected, either geometrically or by linking!"
                          << std::endl;
        }
        if (is_self_link) {
            record_stream << "\n(!!) Found portal crossing inside volume ("
                          << current_rec.volume_idx() << ")!" << std::endl;
        }
        if (is_mixed) {
            record_stream << "\n(!!) Portal crossing involves module surface ("
                          << current_rec.volume_idx() << " <-> "
                          << next_rec.volume_idx() << ")! The second surface in"
                          << " this portal crossing is not a portal!"
                          << std::endl;
        }
        if (is_in_volume and is_valid and not is_self_link and not is_mixed) {

            // Insert into portal trace
            trace_entry lower{current_rec.surface_idx(),
                              current_rec.volume_idx()};
            trace_entry upper{next_rec.surface_idx(), next_rec.volume_idx()};
            portal_trace.emplace_back(lower, upper);
        }
        // Something went wrong
        else {
            // Print search log
            err_stream << "\nError in portal matching:\n" << std::endl;
            err_stream << "volume id\t(intersection info)" << std::endl;

            err_stream << record_stream.str() << std::endl;

            err_stream << "-----\nINFO: Ray terminated at portal x-ing "
                       << (rec + 1) / 2 << ":\n"
                       << current_rec.inters() << " <-> " << next_rec.inters()
                       << std::endl;

            const record rec_front{intersection_records.front()};
            const record rec_back{intersection_records.back()};
            err_stream << "Start volume : " << start_volume << std::endl;
            err_stream << "- first recorded intersection: (sf id:"
                       << rec_front.surface_idx()
                       << ", dist:" << rec_front.dist() << ")," << std::endl;
            err_stream << "- last recorded intersection:  (sf id:"
                       << rec_back.surface_idx() << ", dist:" << rec_back.dist()
                       << "),";

            print_err(err_stream);

            return std::make_tuple(portal_trace, module_trace, error_code);
        }

        // Advance to inspect next pair
        rec += 2u;
    }

    // Look at the last entry, which is a single portal
    const record rec_back{intersection_records.back()};
    if (not rec_back.is_portal()) {
        err_stream << "We don't leave the detector by portal!";
        print_err(err_stream);
    } else {
        trace_entry lower(rec_back.surface_idx(), rec_back.volume_idx());
        trace_entry upper(rec_back.surface_idx(), rec_back.volume_link());
        portal_trace.emplace_back(lower, upper);
    }

    return std::make_tuple(portal_trace, module_trace, error_code);
}

/// Build an adjacency list from intersection traces.
///
/// @tparam portal_trace_type container of portal link pairs
/// @tparam module_trace_type container of module surface links
///
/// @param portal_trace the portal indices and their volume links (in adjacent
///                     portal pairs)
/// @param module_trace the module indices and their volume links
/// @param obj_hashes record which modules/portals were already added
///
/// @return an adjacency list from the traced ray scan of a given geometry.
template <
    dindex invalid_value = dindex_invalid, typename portal_trace_type,
    typename module_trace_type, typename entry_type = std::pair<dindex, dindex>,
    std::enable_if_t<std::is_same_v<typename portal_trace_type::value_type,
                                    std::pair<entry_type, entry_type>>,
                     bool> = true,
    std::enable_if_t<
        std::is_same_v<typename module_trace_type::value_type, entry_type>,
        bool> = true>
inline auto build_adjacency(
    const portal_trace_type &portal_trace,
    const module_trace_type &module_trace,
    std::map<dindex, std::map<dindex, dindex>> &adj_list,
    std::unordered_set<dindex> &obj_hashes) {

    // Every module that was recorded adds a link to the mother volume
    for (const auto &record : module_trace) {
        const auto sf_index = std::get<0>(record);
        const auto vol_index = std::get<1>(record);
        // Check whether we have seen this module in this volume before
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
        if (vol_index_2 != invalid_value) {
            if (obj_hashes.find(pt_index_2) == obj_hashes.end()) {
                adj_list[vol_index_2][vol_index_1]++;
                obj_hashes.insert(pt_index_2);
            }
        }
    }

    return adj_list;
}

/// Build an adjacency list from intersection traces.
///
/// @tparam portal_trace_type container of portal link pairs
/// @tparam module_trace_type container of module links
///
/// @param portal_trace the portal indices and their volume links (in adjacent
///                     portal pairs)
/// @param module_trace the module indices and their volume links
/// @param obj_hashes record which modules/portals were already added
///
/// @return an adjacency list from the traced ray scan of a given geometry.
template <
    dindex invalid_value = dindex_invalid, typename portal_trace_type,
    typename module_trace_type, typename entry_type = std::pair<dindex, dindex>,
    std::enable_if_t<std::is_same_v<typename portal_trace_type::value_type,
                                    std::pair<entry_type, entry_type>>,
                     bool> = true,
    std::enable_if_t<
        std::is_same_v<typename module_trace_type::value_type, entry_type>,
        bool> = true>
inline auto build_adjacency(const portal_trace_type &portal_trace,
                            const module_trace_type &module_trace,
                            dvector<dindex> &adj_matrix,
                            std::unordered_set<dindex> &obj_hashes) {

    const dindex dim = static_cast<dindex>(math::sqrt(adj_matrix.size()));

    // Every module that was recorded adds a link to the mother volume
    for (const auto &record : module_trace) {
        const auto sf_index = std::get<0>(record);
        const auto vol_index = std::get<1>(record);
        // Check whether we have seen this module in this volume before
        if (obj_hashes.find(sf_index) == obj_hashes.end()) {
            adj_matrix[dim * vol_index + vol_index]++;
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
            dindex mat_elem_vol1;
            // Assume the return link for now (filtering out portals that leave
            // world)
            if (vol_index_2 != invalid_value) {
                mat_elem_vol1 = dim * vol_index_1 + vol_index_2;

                if (obj_hashes.find(pt_index_2) == obj_hashes.end()) {
                    adj_matrix[dim * vol_index_2 + vol_index_1]++;
                    obj_hashes.insert(pt_index_2);
                }
            } else {
                mat_elem_vol1 = dim * vol_index_1 + dim - 1;
            }
            adj_matrix[mat_elem_vol1]++;
            obj_hashes.insert(pt_index_1);
        }
    }

    return adj_matrix;
}

}  // namespace detray::detector_scanner
