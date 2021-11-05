/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/utils/ray_gun.hpp"
#include "tests/common/read_geometry.hpp"
#include "tests/common/test_ray_scan.hpp"

using namespace detray;

vecmem::host_memory_resource host_mr;
auto [d, name_map] = read_from_csv(tml_files, host_mr);

namespace __plugin {

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, ray_scan) {

    unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;

    const point3 ori{0., 0., 0.};
    dindex start_index = d.volume_by_pos(ori).index();

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.05 + itheta * (M_PI - 0.1) / theta_steps;
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
            auto [portal_trace, surface_trace] =
                trace_volumes(volume_record, start_index);
            // const auto portal_trace = trace_volumes(volume_record,
            // start_index);

            // All edges made it through the checking
            ASSERT_TRUE(check_connectivity(portal_trace));
        }
    }
}

}  // namespace __plugin

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
