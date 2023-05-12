/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/examples/types.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "tests/common/tools/hash_tree.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

constexpr std::size_t root_hash = 3244;

///
int main() {

    // Toy detector configuration
    vecmem::host_memory_resource host_mr;
    auto det = detray::create_toy_geometry(host_mr);

    // Build the graph
    detray::volume_graph graph(det);
    const auto &adj_mat = graph.adjacency_matrix();

    std::cout << graph.to_string() << std::endl;

    auto geo_checker = detray::hash_tree(adj_mat);

    if (geo_checker.root() == root_hash) {
        std::cout << "Geometry links consistent" << std::endl;
    }
}
