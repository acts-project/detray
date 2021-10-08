/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <utility>

#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {
/**
 * @brief Uses the geometry implementations to walk through their graph-like
 * structure breadth first.
 *
 * This class provides a graph algorithm that walk along the volumes of a given
 * geomtery and uses the portals to check reachability between the volumes.
 *
 * @tparam geometry the type of geometry we want to walk along.
 *
 * @note The geometry has to expose the volume/portal interface.
 */
template <typename geometry, template <typename> class vector_type = dvector>
class geometry_graph {

    public:
    // Objects ids of the geometry
    using object_id = typename geometry::known_objects;

    // Graph nodes
    using node_t = typename geometry::volume_type;

    // Graph edges
    using edge_t = typename geometry::portal;

    /** Default constructor */
    geometry_graph() = delete;

    /** Build from existing nodes and edges, which are provide by the geometry.
     *
     * @param volumes geometry volumes that become the graph nodes
     * @param portals geometry portals link volumes and become edges
     */
    geometry_graph(const vector_type<typename geometry::volume_type> &volumes,
                   const vector_type<typename geometry::portal> &portals)
        : _nodes(volumes), _edges(portals) {
        build();
    }

    /** Default destructor: we don't own anything */
    ~geometry_graph() = default;

    /** @return total number of volumes */
    const size_t n_nodes() const { return _nodes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto nodes() const { return _nodes; }

    /** @return all surfaces/portals in the geometry */
    size_t n_edges() const { return _edges.size(); }

    /** @return all surfaces/portals in the geometry */
    const auto edges() const { return _edges; }

    /**
     * Print geometry if an external name map is provided for the volumes.
     *
     * @param names  Lookup for the names by volume index.
     *
     * @returns the geometry description as a string
     */
    inline const std::string to_string() const {
        std::stringstream ss;
        for (const auto &n : _nodes) {
            ss << "[>>] Node with index " << n.index() << std::endl;
            ss << " -> neighbors: " << n.index() << std::endl;
            const auto &neighbors = adjaceny_list.at(n.index());
            for (const auto &nbr : neighbors) {
                ss << "    -> " << nbr << std::endl;
            }
        }
        return ss.str();
    };

    private:
    /** Walk through the nodes and fill adjacency list. Root node is always at
     * zero.
     */
    void build() {
        for (const auto &n : _nodes) {
            const auto &edge_range = n.template range<object_id::e_portal>();
            vector_type<dindex> neighbors = {};
            std::cout << "on node " << n.index() << std::endl;
            for (size_t edi = edge_range[0]; edi < edge_range[1]; edi++) {
                std::cout << "portal " << edi << std::endl;
                neighbors.push_back(_edges[edi].volume());
            }
            std::cout << "found neighbors " << neighbors.size() << std::endl;
            adjaceny_list[n.index()] = neighbors;
            std::cout << "adj list size: " << adjaceny_list.size() << std::endl;
        }
    }

    void bfs() {
        // Start at root
        const auto &current = _nodes.front();

        // Nodes to be visited
        std::queue<node_t> node_queue;
        node_queue.push(current);

        // Visit adjacent nodes
        do {
            const auto &neighbors = adjaceny_list.at(current.index());
            current = node_queue.pop();
        } while (not node_queue.empty());
    }

    /** Graph nodes */
    const vector_type<node_t> &_nodes;

    /** Graph edges */
    const vector_type<edge_t> &_edges;

    /** Volume indices, where reachability was established. The source link
     * is encoded as the position of the neigbor vector in the outer vector.
     */
    std::map<dindex, vector_type<dindex>> adjaceny_list = {};
};

}  // namespace detray