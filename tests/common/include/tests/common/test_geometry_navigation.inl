/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "core/detector.hpp"
#include "io/csv_io.hpp"
#include "utils/ray_gun.hpp"

using namespace detray;

/** Read the detector from file */
auto read_detector() {
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string name = "tml";
    std::string surface_file = data_directory + "tml.csv";
    std::string layer_volume_file = data_directory + "tml-layer-volumes.csv";
    std::string surface_grid_file = data_directory + "tml-surface-grids.csv";
    std::string surface_grid_entries_file = "";

    /*std::string name = "odd";
    std::string surface_file = data_directory + std::string("odd.csv");
    std::string layer_volume_file =
        data_directory + std::string("odd-layer-volumes.csv");
    std::string surface_grid_file =
        data_directory + std::string("odd-surface-grids.csv");
    std::string surface_grid_entries_file = "";*/

    return detector_from_csv<>(name, surface_file, layer_volume_file,
                               surface_grid_file, surface_grid_entries_file);
};

/** Check if a set of volume index pairs form a trace.
 *
 * @param volume_records the recorded portal crossings between volumes
 * @param start_volume where the ray started
 *
 * @note Empty traces pass automatically.
 *
 * @return true if the volumes form a connected chain.
 */
inline bool trace_volumes(std::set<std::pair<dindex, dindex>> volume_records,
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
inline auto check_connectivity(const record_type &volume_record,
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

auto d = read_detector();

namespace __plugin {

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, ray_scan) {

    unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;
    const unsigned int itest = 10000;

    const point3 ori{0., 0., 0.};
    dindex start_index = d.volume_by_pos(ori).index();

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const point3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                             cos_theta};

            const auto volume_record = shoot_ray(d, ori, dir);
            const auto volume_connections =
                check_connectivity(volume_record, start_index);

            // All edges made it through the checking
            ASSERT_TRUE(trace_volumes(volume_connections));
        }
    }
}

}  // namespace __plugin

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
