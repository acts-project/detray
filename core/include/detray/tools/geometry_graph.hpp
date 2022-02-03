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

#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** Placeholder struct for the implementation of an inspection function
 * that will be executed, when a node is visited.
 */
template <typename node_t>
struct void_node_inspector {
    bool operator()(const node_t &n) { return true; }
};

/** Placeholder struct for an action while walking through the graph.*/
template <typename node_t>
struct void_actor {
    void operator()(const node_t & /*n*/, const dindex_range & /*edge_range*/) {
    }
};

/**
 * @brief Uses the geometry implementations to walk through their graph-like
 * structure breadth first.
 *
 * This class provides a graph algorithm that walk along the volumes of a given
 * geomtery and uses the portals to check reachability between the volumes.
 *
 * @tparam geometry the type of geometry we want to walk along.
 * @tparam node_inspector the type of inspection to perform when a node is
 *         visited
 *
 * @note The geometry has to expose the volume/portal interface.
 */
template <typename geometry,
          typename node_inspector =
              void_node_inspector<typename geometry::volume_type>,
          template <typename...> class vector_t = dvector>
class geometry_graph {

    public:
    // Objects ids of the geometry
    using object_id = typename geometry::object_id;

    // Graph nodes
    using node_t = typename geometry::volume_type;

    // Graph edges
    using edge_t = typename geometry::surface_type;

    /** Default constructor */
    geometry_graph() = delete;

    /** Build from existing nodes and edges, which are provide by the geometry.
     *
     * @param volumes geometry volumes that become the graph nodes
     * @param portals geometry portals link volumes and become edges
     */
    geometry_graph(const vector_t<typename geometry::volume_type> &volumes,
                   const vector_t<typename geometry::surface_type> &portals)
        : _nodes(volumes), _edges(portals) {
        build();
    }

    /** Default destructor: we don't own anything */
    ~geometry_graph() = default;

    /** @return number of volumes */
    size_t n_nodes() const { return _nodes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto &nodes() const { return _nodes; }

    /** @return number of surfaces/portals in the geometry */
    size_t n_edges() const { return _edges.size(); }

    /** @return all surfaces/portals in the geometry */
    const auto &edges() const { return _edges; }

    /** @return graph adjacency */
    auto &adjacency_list() const { return adj_list; }

    /** Walks breadth first through the geometry objects. */
    template <typename action_t = void_actor<node_t>>
    void bfs(action_t actor = {}) const {
        // Do node inspection
        const auto inspector = node_inspector();

        node_t const *current = nullptr;
        vector_t<bool> visited(_nodes.size(), false);

        // Nodes to be visited. Start at first node
        std::queue<node_t const *> node_queue;
        node_queue.push(&(_nodes[0]));

        // Visit adjacent nodes and check current one
        while (not node_queue.empty()) {
            // Inspect
            current = node_queue.front();
            if (visited[current->index()]) {
                node_queue.pop();
                continue;
            }
            inspector(*current);

            // Visit neighbors and perform action
            actor(*current, current->template range<object_id::e_portal>());

            // Add neightbors to queue
            for (const auto &edg : range(_edges, *current)) {
                // Retrieve the node index the edge points to
                dindex nbr = std::get<0>(edg.edge());
                // If not leaving world and if not visited, enqueue the node
                if ((nbr != dindex_invalid and nbr > 0) and not visited[nbr]) {
                    node_queue.push(&(_nodes[nbr]));
                }
            }

            // This node has been visited
            visited[current->index()] = true;
            node_queue.pop();
        }
    }

    /** @returns the linking description as a string */
    inline const std::string to_string() const {
        std::stringstream ss;
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
        for (const auto &n : _nodes) {
            ss << "[>>] Node with index " << n.index() << std::endl;
            ss << " -> neighbors: " << std::endl;
            const auto &neighbors = adj_list.at(n.index());
            for (const auto &nbr : neighbors) {
                ss << "    -> " << print_neighbor(nbr) << std::endl;
            }
        }
        return ss.str();
    };

    private:
    /** Go through the nodes and fill adjacency list. Root node is always at
     * zero.
     */
    void build() {
        for (const auto &n : _nodes) {
            // Count the number of edges for a particluar neighbor
            std::map<dindex, dindex> neighbors = {};

            // Only works for non batched geometries
            for (const auto &edg : range(_edges, n)) {
                neighbors[std::get<0>(edg.edge())]++;
            }
            adj_list[n.index()] = neighbors;
        }
    }

    /** Graph nodes */
    const vector_t<node_t> &_nodes;

    /** Graph edges */
    const vector_t<edge_t> &_edges;

    /**
     *  The index of the nodes neighbors and a count of edges is kept in the
     *  inner map.
     */
    std::map<dindex, std::map<dindex, dindex>> adj_list = {};
};

}  // namespace detray