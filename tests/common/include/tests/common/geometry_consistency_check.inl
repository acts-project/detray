/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>

#include "tests/common/read_geometry.hpp"
#include "tools/geometry_graph.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

// Adhoc geometry type. Or use any other geometry
template <typename volume_t, typename object_t,
          template <typename> class vector_type = dvector>
struct dummy_geometry {
    using volume_type = volume_t;
    using portal = object_t;

    // Known primitives
    enum known_objects : bool {
        e_surface = true,
        e_portal = false,
        e_any = false,  // defaults to portal
    };

    dummy_geometry(const vector_type<volume_t> &volumes,
                   const vector_type<object_t> &objects)
        : _volumes(volumes), _objects(objects) {}

    const vector_type<volume_t> &_volumes;
    const vector_type<object_t> &_objects;
};

// Tests the consistency of the toy geometry implementation. In principle,
// every geometry can be checked this way.
TEST(ALGEBRA_PLUGIN, geometry_consistency) {

    // Build the geometry (modeled as a unified index geometry)
    auto [volumes, surfaces, transforms, cylinders, discs, rectangles] =
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

    std::cout << g.to_string() << std::endl;
    // std::cout << "Walking through geometry: " << std::endl;
    g.bfs();

    const auto &adj = g.adjacency_list();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
