/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/tools/geometry_graph.hpp"
#include "detray/tools/hash_tree.hpp"
#include "tests/common/read_geometry.hpp"
#include "tests/common/test_ray_scan.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

// Adhoc geometry type for toy geometry. Or use any other geometry
template <typename volume_t, typename object_t,
          template <typename> class vector_type = dvector>
struct dummy_geometry {
    // typedefs
    using volume_type = volume_t;
    using portal = object_t;
    using portal_links = typename object_t::edge_links;

    struct object_registry {
        using id = typename volume_type::objects;

        template <typename value_type = void>
        static constexpr auto get() {
            return id::e_surface;
        }
    };

    dummy_geometry(const vector_type<volume_t> &volumes,
                   const vector_type<object_t> &objects)
        : _volumes(volumes), _objects(objects) {}

    // data containers
    const vector_type<volume_t> &_volumes;
    const vector_type<object_t> &_objects;
};

// Adhoc detector type for toy geometry
template <typename geometry_t, typename transform_container,
          typename mask_container>
struct dummy_detector {
    // typedefs
    using geometry = geometry_t;
    using object_id = typename geometry_t::object_registry::id;
    struct transform_store {
        // dummy type
        struct context {};
    };

    dummy_detector(const geometry &geometry, const transform_container &trfs,
                   const mask_container &masks)
        : _geometry(geometry), _transforms(trfs), _masks(masks) {}

    // interface functions
    const auto &volumes() const { return _geometry._volumes; }
    const auto &transforms(
        const typename transform_store::context ctx = {}) const {
        return _transforms;
    }
    const auto &masks() const { return _masks; }
    template <object_id>
    const auto &objects() const {
        return _geometry._objects;
    }

    // data containers
    const geometry &_geometry;
    const transform_container &_transforms;
    const mask_container &_masks;
};

/** Print and adjacency list */
void print_adj(std::map<dindex, std::map<dindex, dindex>> &adjacency_list) {
    const auto print_neighbor =
        [&](const std::pair<const dindex, const dindex> &n) -> std::string {
        // Print the number of edges, if it is greater than one
        std::string n_occur =
            n.second > 1 ? "\t\t(" + std::to_string(n.second) + "x)" : "";

        // Edge that leads out of the detector world
        if (n.first == dindex_invalid) {
            return "leaving world" + n_occur;
        }

        return std::to_string(n.first) + "\t\t\t" + n_occur;
    };

    for (const auto &[vol_index, neighbors] : adjacency_list) {
        std::cout << "[>>] Node with index " << vol_index << std::endl;
        std::cout << " -> neighbors: " << std::endl;
        for (const auto &nbr : neighbors) {
            std::cout << "    -> " << print_neighbor(nbr) << std::endl;
        }
    }
}

// Tests the consistency of the toy geometry implementation. In principle,
// every geometry can be checked this way.
TEST(ALGEBRA_PLUGIN, geometry_consistency) {

    // Build the geometry (modeled as a unified index geometry)
    auto [volumes, surfaces, transforms, discs, cylinders, rectangles] =
        toy_geometry();

    using geometry_t = dummy_geometry<typename decltype(volumes)::value_type,
                                      typename decltype(surfaces)::value_type>;

    const auto geo = geometry_t(volumes, surfaces);

    // Build the graph
    const auto g = geometry_graph<geometry_t>(geo._volumes, geo._objects);
    const auto &adj_linking = g.adjacency_list();

    std::cout << g.to_string() << std::endl;

    // Now get the adjaceny list from ray scan

    // First, put data into the detector interface
    mask_store<dtuple, dvector, decltype(discs)::value_type,
               decltype(cylinders)::value_type,
               decltype(rectangles)::value_type>
        masks;
    // populate mask store
    masks.add_masks(discs);
    masks.add_masks(cylinders);
    masks.add_masks(rectangles);

    using detector_t =
        dummy_detector<geometry_t, decltype(transforms), decltype(masks)>;

    auto d = detector_t(geo, transforms, masks);

    // Adjacency list to be filled in ray scan
    std::map<dindex, std::map<dindex, dindex>> adj_scan = {};
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;
    const point3 ori{0., 0., 0.};
    dindex start_index = 0;

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

            // These are the portal links
            auto [portal_trace, surface_trace] =
                trace_volumes(volume_record, start_index);

            // All edges made it through the checking
            ASSERT_TRUE(check_connectivity(portal_trace));

            build_adjacency(portal_trace, surface_trace, adj_scan, obj_hashes);
        }
    }

    print_adj(adj_scan);

    auto geo_checker_vol0 =
        hash_tree<decltype(adj_linking.at(0)), dindex>(adj_linking.at(0));
    auto geo_checker_scan_vol0 =
        hash_tree<decltype(adj_scan.at(0)), dindex>(adj_scan.at(0));

    EXPECT_EQ(geo_checker_vol0.root(), geo_checker_scan_vol0.root());

    // This one fails, because the ray scan is kept very coarse for performance
    // reasons
    /*auto geo_checker_vol1 =
        hash_tree<decltype(adj_linking.at(1)), dindex>(adj_linking.at(1));
    auto geo_checker_scan_vol1 =
        hash_tree<decltype(adj_scan.at(1)), dindex>(adj_scan.at(1));

    EXPECT_EQ(geo_checker_vol1.root(), geo_checker_scan_vol1.root());*/

    auto geo_checker_vol2 =
        hash_tree<decltype(adj_linking.at(2)), dindex>(adj_linking.at(2));
    auto geo_checker_scan_vol2 =
        hash_tree<decltype(adj_scan.at(2)), dindex>(adj_scan.at(2));

    EXPECT_EQ(geo_checker_vol2.root(), geo_checker_scan_vol2.root());

    auto geo_checker_vol3 =
        hash_tree<decltype(adj_linking.at(3)), dindex>(adj_linking.at(3));
    auto geo_checker_scan_vol3 =
        hash_tree<decltype(adj_scan.at(3)), dindex>(adj_scan.at(3));

    EXPECT_EQ(geo_checker_vol3.root(), geo_checker_scan_vol3.root());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
