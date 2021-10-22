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
    void operator()(const node_t &n, const dindex_range &edge_range) {}
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
          template <typename> class vector_type = dvector>
class geometry_graph {

    public:
    // Objects ids of the geometry
    using object_id = typename geometry::object_id;

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

    /** @return number of volumes */
    const size_t n_nodes() const { return _nodes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto nodes() const { return _nodes; }

    /** @return number of surfaces/portals in the geometry */
    size_t n_edges() const { return _edges.size(); }

    /** @return all surfaces/portals in the geometry */
    const auto edges() const { return _edges; }

    /** @return graph adjacency */
    const auto adjacency_list() const { return adjaceny_list; }

    /** Walks breadth first through the geometry objects. */
    template <typename action_t = void_actor<node_t>>
    void bfs(action_t actor = {}) const {
        // Do node inspection
        const auto inspector = node_inspector();

        node_t const *current = nullptr;
        vector_type<bool> visited(_nodes.size(), false);

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
            const auto &edge_range =
                current->template range<object_id::e_portal>();
            actor(*current, edge_range);

            // Add neightbors to queue
            for (size_t edi = edge_range[0]; edi < edge_range[1]; edi++) {
                // Retrieve the node index the edge points to
                dindex nbr = std::get<0>(_edges[edi].edge());
                // If not visited, enqueue the node
                if (not visited[nbr]) {
                    node_queue.push(&(_nodes[nbr]));
                }
            }

            // This volume has been visited
            visited[current->index()] = true;
            node_queue.pop();
        }
    }

    /** @returns the linking description as a string */
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
    /** Go through the nodes and fill adjacency list. Root node is always at
     * zero.
     */
    void build() {
        for (const auto &n : _nodes) {
            const auto &edge_range = n.template range<object_id::e_portal>();
            vector_type<dindex> neighbors = {};

            for (size_t edi = edge_range[0]; edi < edge_range[1]; edi++) {
                neighbors.push_back(std::get<0>(_edges[edi].edge()));
            }
            adjaceny_list[n.index()] = neighbors;
        }
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