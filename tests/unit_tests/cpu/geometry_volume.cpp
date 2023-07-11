/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/detail/volume_descriptor.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/test/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// TODO: Move these into the test defs
namespace {

constexpr detray::scalar tol{1e-6f};

// geo object ids for testing
enum geo_objects : unsigned int {
    e_sensitive = 0,
    e_portal = 1,
    e_size = 2,
    e_all = e_size,
};

// surface finder ids for testing
enum sf_finder_ids : unsigned int {
    e_default = 0,
    e_grid = 1,
};

}  // namespace

// This tests the detector volume descriptor and its many links
GTEST_TEST(detray_geometry, detector_volume_descriptor) {
    using namespace detray;

    using sf_finder_link_t = dtyped_index<sf_finder_ids, dindex>;
    using volume_t = volume_descriptor<geo_objects, sf_finder_link_t>;

    // Check construction, setters and getters
    volume_t v1(volume_id::e_cylinder);
    v1.set_index(12345u);
    v1.template set_link<geo_objects::e_portal>({sf_finder_ids::e_default, 1u});
    v1.template set_link<geo_objects::e_sensitive>(
        {sf_finder_ids::e_grid, 12u});

    ASSERT_TRUE(v1.id() == volume_id::e_cylinder);
    ASSERT_TRUE(v1.index() == 12345u);
    ASSERT_TRUE(v1.template link<geo_objects::e_portal>().id() ==
                sf_finder_ids::e_default);
    ASSERT_TRUE(v1.template link<geo_objects::e_portal>().index() == 1u);
    ASSERT_TRUE(v1.template link<geo_objects::e_sensitive>().id() ==
                sf_finder_ids::e_grid);
    ASSERT_TRUE(v1.template link<geo_objects::e_sensitive>().index() == 12u);

    // Check copy constructor
    const auto v2 = volume_t(v1);
    ASSERT_EQ(v2.id(), volume_id::e_cylinder);
    ASSERT_EQ(v2.index(), 12345u);
    ASSERT_TRUE(v2.template link<geo_objects::e_portal>().id() ==
                sf_finder_ids::e_default);
    ASSERT_TRUE(v2.template link<geo_objects::e_portal>().index() == 1u);
    ASSERT_TRUE(v2.template link<geo_objects::e_sensitive>().id() ==
                sf_finder_ids::e_grid);
    ASSERT_TRUE(v2.template link<geo_objects::e_sensitive>().index() == 12u);
}

// This tests the detector volume interface class
GTEST_TEST(detray_geometry, detector_volume) {
    using namespace detray;

    using detector_t = detector<toy_metadata<>>;

    vecmem::host_memory_resource host_mr;
    const auto toy_det = create_toy_geometry(host_mr);

    // This is the beampipe volume
    const auto vol = toy_det.volume_by_index(0u);

    // TODO: Add beampipe name once name map PR is merged
    ASSERT_EQ(vol.id(), volume_id::e_cylinder);
    ASSERT_EQ(vol.index(), 0u);
    // Identity
    ASSERT_TRUE(vol.transform() == detector_t::transform3{});
    const auto c = vol.center();
    ASSERT_NEAR(c[0], 0.f, tol);
    ASSERT_NEAR(c[1], 0.f, tol);
    ASSERT_NEAR(c[2], 0.f, tol);
    // The beampipe + all portals
    ASSERT_EQ(vol.n_max_candidates(), 16u);

    // Check the portal access
    const auto portals = vol.portals();
    ASSERT_EQ(portals.size(), 16u);
    for (const auto [idx, pt_desc] : detray::views::enumerate(portals)) {
        const auto portal = surface{toy_det, pt_desc};

        // In this volume, also the beampipe will be among the portals
        ASSERT_TRUE(portal.is_portal() or portal.is_passive())
            << "Not a portal/passive: " << portal;
        ASSERT_EQ(portal.volume(), 0u) << "Wrong volume: " << portal;
        ASSERT_EQ(portal.index(), idx) << "Wrong index: " << portal;
    }
}
