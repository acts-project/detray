/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tools/intersection_kernel.hpp"

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

/** Intersect all portals in a detector with a given ray.
 *
 * @tparam detector_type the type of the detector
 *
 * @param d the detector
 * @param origin the origin of the ray in global coordinates
 * @param direction the direction of the ray in global coordinater
 *
 * @return a sorted vector of volume indices with the corresponding
 *         intersections of the portals that were encountered
 */
template <typename detector_type>
auto shoot_ray(const detector_type &d, const std::pair<point3, point3> origin,
               const std::pair<point3, point3> direction) {

    using object_id = typename detector_type::objects;
    using portal_links = typename detector_type::geometry::portal_links;
    using detray_context = typename detector_type::transform_store::context;

    detray_context default_context;

    track<detray_context> ray;
    ray.pos = origin.first;
    ray.dir = direction.first;

    std::vector<std::pair<dindex, intersection>> volume_record;

    const auto &transforms = d.transforms(default_context);
    const auto &portals = d.portals();
    const auto &masks = d.masks();
    // Loop over volumes
    for (const auto &v : d.volumes()) {
        // Record the portals the ray intersects
        const auto &pt_range = v.template range<object_id::e_portal>();
        portal_links links = {};
        for (size_t pi = pt_range[0]; pi < pt_range[1]; pi++) {
            auto pti = intersect(ray, portals[pi], transforms, masks, links);

            if (pti.status == intersection_status::e_inside &&
                pti.direction == intersection_direction::e_along) {
                volume_record.emplace_back(v.index(), pti);
            }
        }
    }

    // Sort by distance to origin of the ray and then volume index
    auto sort_path = [&](std::pair<dindex, intersection> a,
                         std::pair<dindex, intersection> b) -> bool {
        return (a.second.path == b.second.path) ? (a.first < b.first)
                                                : (a.second < b.second);
    };
    std::sort(volume_record.begin(), volume_record.end(), sort_path);

    return volume_record;
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
auto check_record(const record_type &volume_record, dindex start_volume = 0) {
    std::set<std::pair<dindex, dindex>> valid_volumes = {};

    // Always read 2 elements from the sorted records vector
    struct record_doublet {
        // Two volume index, portal intersections
        const typename record_type::value_type &lower, upper;

        // Is this doublet connected via the portal intersection?
        inline bool operator()() const { return lower.second == upper.second; }
    };

    // Don't look a start and end volume, as they are not connected at one side
    for (size_t rec = 0; rec < volume_record.size() - 1; rec += 2) {
        // Get 2 possibly connected entries
        record_doublet doublet = {.lower = volume_record[rec],
                                  .upper = volume_record[rec + 1]};

        /*std::cout << doublet.lower.first << " (" << doublet.lower.second.path
        << ")" << std::endl; std::cout << doublet.upper.first << " (" <<
        doublet.upper.second.path << ")" << std::endl; std::cout << std::endl;*/

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
                check_record(volume_record, start_index);

            // All edges made it through the checking
            check_result =
                (volume_record.size() - 1) / 2 == volume_connections.size();
        }
    }
    // Did all rays pass?
    ASSERT_TRUE(check_result);
}

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, random_rays) {}

}  // namespace __plugin

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
