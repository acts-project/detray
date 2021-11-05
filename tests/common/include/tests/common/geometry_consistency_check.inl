/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>

#include "tests/common/read_geometry.hpp"
#include "tests/common/test_ray_scan.hpp"
#include "detray/tools/geometry_graph.hpp"
//#include <vecmem/memory/host_memory_resource.hpp>

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

// Tests the consistency of the toy geometry implementation. In principle,
// every geometry can be checked this way.
TEST(ALGEBRA_PLUGIN, geometry_consistency) {

    // Build the geometry (modeled as a unified index geometry)
    auto [volumes, surfaces, transforms, discs, cylinders, rectangles] =
        toy_geometry();

    using geometry_t = dummy_geometry<typename decltype(volumes)::value_type,
                                      typename decltype(surfaces)::value_type>;

    const auto geo = geometry_t(volumes, surfaces);

    /// Prints linking information for every node when visited
    struct volume_printout {
        void operator()(const geometry_t::volume_type &n) const {
            std::cout << "On volume: " << n.index() << std::endl;
        }
    };

    using graph = geometry_graph<geometry_t, volume_printout>;

    // Build the graph
    graph g = graph(geo._volumes, geo._objects);
    const auto &adj = g.adjacency_list();

    std::cout << g.to_string() << std::endl;

    /*vecmem::host_memory_resource host_mr;
    auto [d, name_map] = read_from_csv(odd_files, host_mr);
    auto odd_graph = geometry_graph<decltype(d)::geometry,
    volume_printout>(d.volumes(), d.template
    objects<decltype(d)::geometry::object_registry::id::e_portal>());

    std::cout << odd_graph.to_string() << std::endl;*/

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
            const auto volume_trace = trace_volumes(volume_record, start_index);

            // All edges made it through the checking
            ASSERT_TRUE(check_connectivity(volume_trace));
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
