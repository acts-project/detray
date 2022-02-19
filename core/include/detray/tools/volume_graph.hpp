/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cstddef>
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
 * @tparam detector_t the type of geometry we want to walk along.
 * @tparam node_inspector the type of inspection to perform when a node is
 *         visited
 *
 * @note The detector has to expose the volume/portal interface.
 */
template <typename detector_t,
          typename node_inspector =
              void_node_inspector<typename detector_t::volume_type>,
          template <typename...> class vector_t = dvector>
class volume_graph {

    public:
    // Objects ids of the geometry
    using object_defs = typename detector_t::objects;
    using volume_container_t = vector_t<typename detector_t::volume_type>;
    using surface_container_t = vector_t<typename detector_t::surface_type>;
    using mask_container_t = typename detector_t::mask_container;

    /** Builds a graph node from the detector collections on the fly. */
    struct node_collection {
        struct node {

            node(const typename detector_t::volume_type &volume,
                 const surface_container_t &surfaces)
                : _idx(volume.index()) {
                for (const auto &sf : range(surfaces, volume)) {
                    _edge_links.push_back(sf.mask());
                }
            }

            dindex index() const { return _idx; }
            const auto &edges() const { return _edge_links; }

            dindex _idx;
            vector_t<typename detector_t::surface_type::mask_link> _edge_links;
        };

        struct iterator {
            using value_type = node;
            using volume_iter =
                decltype(std::begin(std::declval<volume_container_t>()));

            iterator(volume_iter &&vol_itr, const surface_container_t &surfaces)
                : _vol_itr(vol_itr), _surfaces(surfaces) {}

            node operator*() { return node(*_vol_itr, _surfaces); }

            // Prefix increment
            iterator &operator++() {
                ++_vol_itr;
                return *this;
            }

            bool operator==(const iterator &rhs) {
                return _vol_itr == rhs._vol_itr and
                       &_surfaces == &rhs._surfaces;
            };

            bool operator!=(const iterator &rhs) {
                return _vol_itr != rhs._vol_itr or &_surfaces != &rhs._surfaces;
            };

            volume_iter _vol_itr;
            const surface_container_t &_surfaces;
        };

        using value_type = node;

        node_collection(const volume_container_t &volumes,
                        const surface_container_t &surfaces)
            : _volumes(volumes), _surfaces(surfaces) {}

        std::size_t size() const { return _volumes.size(); }
        node front() const { return node(_volumes.front(), _surfaces); }
        iterator begin() const { return iterator(_volumes.begin(), _surfaces); }
        iterator end() const { return iterator(_volumes.end(), _surfaces); }

        const volume_container_t &_volumes;
        const surface_container_t &_surfaces;
    };

    /** Builds a graph edge from the detector mask collection on the fly. */
    struct edge_collection {
        using mask_edge_t = typename detector_t::surface_type::edge_type;
        using mask_link_t = typename detector_t::surface_type::mask_link;

        struct edge {
            edge(const dindex volume_id, const dindex link)
                : _from(volume_id), _to(link) {}

            dindex from() const { return _from; }
            dindex to() const { return _to; }

            dindex _from, _to;
        };

        struct iterator {
            using value_type = edge;
            using edge_iter =
                decltype(std::begin(std::declval<vector_t<edge>>()));

            iterator(const dindex volume_id, const mask_link_t &mask_link,
                     const edge_collection &edges) {
                build_edges_vector(volume_id, mask_link, edges.get_container());
                _itr = _edges.begin();
            }

            edge operator*() { return *_itr; }

            // Prefix increment
            iterator &operator++() {
                ++_itr;
                return *this;
            }

            iterator &begin() {
                _itr = _edges.begin();
                return *this;
            }
            iterator &end() {
                _itr = _edges.end();
                return *this;
            }

            bool operator==(const iterator &rhs) {
                return _itr == rhs._itr and &_edges == &rhs._edges;
            };

            bool operator!=(const iterator &rhs) {
                return _itr != rhs._itr and &_edges != &rhs._edges;
            };

            template <std::size_t current_id = 0>
            inline void build_edges_vector(const dindex volume_id,
                                           const mask_link_t &mask_link,
                                           const mask_container_t &masks) {

                if (detail::get<0>(mask_link) == current_id) {
                    // Get the mask group that will be updated
                    const auto &mask_group = masks.template group<current_id>();
                    const auto mask_range = detail::get<1>(mask_link);
                    for (const auto &mask : range(mask_group, mask_range)) {
                        _edges.emplace_back(volume_id, mask.volume_link());
                    }
                }

                // Next mask type
                using mask_defs = typename detector_t::surface_type::mask_defs;
                if constexpr (current_id < mask_defs::n_types - 1) {
                    return build_edges_vector<current_id + 1>(volume_id,
                                                              mask_link, masks);
                }
            }

            edge_iter _itr;
            vector_t<edge> _edges;
        };

        using value_type = edge;

        edge_collection(const typename detector_t::mask_container &masks)
            : _edges(masks) {}

        std::size_t size() const { return dindex_invalid; }
        const mask_container_t &get_container() const { return _edges; }

        const mask_container_t &_edges;
    };

    // Graph nodes
    using node_type = typename node_collection::value_type;
    using node_iter = typename node_collection::iterator;
    // Graph edges
    using edge_type = typename edge_collection::value_type;
    using edge_iter = typename edge_collection::iterator;

    /** Default constructor */
    volume_graph() = delete;

    /** Build from existing nodes and edges, which are provide by the detector.
     *
     * @param det provides: geometry volumes that become the graph nodes,
     * surfaces which are needed to index the correct masks and the
     * masks that link to volumes and become graph edges
     */
    volume_graph(const detector_t &det)
        : _nodes(det.volumes(), det.surfaces()), _edges(det.mask_store()) {
        build();
    }

    /** Default destructor: we don't own anything */
    ~volume_graph() = default;

    /** @return number of volumes */
    size_t n_nodes() const { return _nodes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto &nodes() const { return _nodes; }

    /** @return number of surfaces/portals in the geometry */
    // size_t n_edges() const { return _edges.size(); }

    /** @return all surfaces/portals in the geometry */
    const auto &edges() const { return _edges; }

    /** @return graph adjacency */
    auto &adjacency_list() const { return adj_list; }

    /** Walks breadth first through the geometry objects. */
    /*template <typename action_t = void_actor<node_type>>
    void bfs(action_t actor = {}) const {
        // Do node inspection
        const auto inspector = node_inspector();

        node_type const *current = nullptr;
        vector_t<bool> visited(_nodes.size(), false);

        // Nodes to be visited. Start at first node
        std::queue<node_type const *> node_queue;
        node first_node = _nodes.front();
        node_queue.push(&(first_node));

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
            actor(*current, current->template range<object_defs::e_portal>());

            // Add neightbors to queue
            for (const auto &edg_link : current->edges()) {
                for (const auto &edg :
    edge_collection::iterator(current->index(), edg_link, _edges)) { dindex nbr
    = edg.to();
                    // If not leaving world and if not visited, enqueue the node
                    if ((nbr != dindex_invalid and nbr > 0) and not
    visited[nbr]) { node_queue.push(&(_nodes[nbr]));
                    }
                }
            }

            // This node has been visited
            visited[current->index()] = true;
            node_queue.pop();
        }
    }*/

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
     *  zero.
     */
    void build() {
        for (const auto &n : _nodes) {
            // Count the number of edges for a particluar neighbor
            std::map<dindex, dindex> neighbors = {};

            for (const auto &edg_link : n.edges()) {
                for (const auto edg : edge_iter(n.index(), edg_link, _edges)) {
                    neighbors[edg.to()]++;
                }
            }
            adj_list[n.index()] = neighbors;
        }
    }

    /** Graph nodes */
    node_collection _nodes;

    /** Graph edges */
    edge_collection _edges;

    /**
     *  The index of the nodes neighbors and a count of edges is kept in the
     *  inner map.
     */
    std::map<dindex, std::map<dindex, dindex>> adj_list = {};
};

}  // namespace detray