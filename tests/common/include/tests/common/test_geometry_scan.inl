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
    std::string surfaces = data_directory + "tml.csv";
    std::string volumes = data_directory + "tml-layer-volumes.csv";
    std::string grids = data_directory + "tml-surface-grids.csv";
    std::string grid_entries = "";
    // std::map<dindex, std::string> name_map{};

    /*return detray::detector_from_csv<>(name, surfaces, volumes, grids,
                                       grid_entries, name_map);*/
    return detray::detector_from_csv<>(name, surfaces, volumes, grids,
                                       grid_entries);
};

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
inline auto check_connectivity(const record_type &volume_record, dindex start_volume = 0) {
    std::set<std::pair<dindex, dindex>> valid_volumes = {};

    // Always read 2 elements from the sorted records vector
    struct record_doublet {
        // Two volume index, portal intersections
        const typename record_type::value_type &lower, upper;

        // Is this doublet connected via the portal intersection?
        inline bool operator()() const { return lower.second == upper.second; }
    };

    if ((volume_record.size() - 1) % 2) {
        return valid_volumes;
    }

    // Don't look a start and end volume, as they are not connected at one side
    for (size_t rec = 1; rec < volume_record.size() - 1; rec += 2) {

        // Get 2 possibly connected entries
        record_doublet doublet = {.lower = volume_record[rec],
                                  .upper = volume_record[rec + 1]};

        if (doublet()) {
            // Insert into set of edges
            valid_volumes.emplace(doublet.lower.first, doublet.upper.first);
        }
        // Something went wrong
        else {
            std::cerr << "<<<<<<<<<<<<<<<" << std::endl;

            std::cerr << "Ray terminated at portal x-ing " << (rec + 1) / 2
                      << ": (" << doublet.lower.first << ", "
                      << doublet.lower.second.path << ") <-> ("
                      << doublet.upper.first << ", "
                      << doublet.upper.second.path << ")" << std::endl;

            std::cerr << "Start volume : " << start_volume
                      << ", leaves world at: (" << volume_record.front().first
                      << ", " << volume_record.front().second.path << "), ("
                      << volume_record.back().first << ", "
                      << volume_record.back().second.path << ")" << std::endl;

            if (volume_record.front().second == doublet.lower.second) {
                std::cerr << "=> Ray didn't leave start volume at "
                          << volume_record.front().first << std::endl;
            }
            if (volume_record.back().second == doublet.upper.second) {
                std::cerr << "=> Ray didn't leave world at "
                          << volume_record.back().first << std::endl;
            }

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

    const auto ori = std::make_pair<point3, point3>({0., 0., 0.}, {0., 0., 0.});
    dindex start_index = d.volume_by_pos(ori.first).index();
    bool check_result = false;

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
            const auto dir = std::make_pair<point3, point3>(
                {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta},
                {0., 0., 0.});

            const auto volume_record = shoot_ray(d, ori, dir);
            const auto volume_connections =
                check_connectivity(volume_record, start_index);

            // All edges made it through the checking
            check_result =
                (volume_record.size() - 1) / 2 == volume_connections.size();
        }
    }
    // Did all rays pass?
    ASSERT_TRUE(check_result);
}

}  // namespace __plugin

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
