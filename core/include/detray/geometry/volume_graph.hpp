/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cstddef>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <utility>

#include "detray/definitions/indexing.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

/// @brief Placeholder struct for the implementation of an inspection function
/// that will be executed, when a node is visited.
template <typename node_t>
struct void_node_inspector {
    bool operator()(const node_t & /*n*/) { return true; }
};

/// @brief Placeholder struct for an action while walking through the graph.
template <typename node_t>
struct void_actor {
    void operator()(const node_t & /*n*/, const dindex_range & /*edge_range*/) {
    }
};

/// @brief Uses the geometry implementations to walk through their graph-like
/// structure breadth first.
///
/// This class provides a graph algorithm that walk along the volumes of a given
/// geomtery and uses the portals to check reachability between the volumes.
///
/// @tparam detector_t the type of geometry we want to walk along.
/// @tparam node_inspector the type of inspection to perform when a node is
///         visited
///
/// @note The detector has to expose the volume/portal interface.
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

    /// @brief Builds a graph node from the detector collections on the fly.
    ///
    /// The node collection builds a graph node from a volume and the mask
    /// links in its surfaces, which are half edges of the graph. Every node
    /// has its index in the detector volume container and a vector of half
    /// edges.
    struct node_collection {

        /// A node in the graph has an index (volume index) and a collction of
        /// edge that belong to it (mask link of every surface in the volume).
        /// One mask link is a half edge in the graph.
        struct node {

            /// Constructor from a detectors volume and surface collections
            node(const typename detector_t::volume_type &volume,
                 const surface_container_t &surfaces)
                : _idx(volume.index()) {
                for (const auto &sf : range(surfaces, volume)) {
                    _half_edges.push_back(sf.mask());
                }
            }

            /// @returns volume index of the node
            dindex index() const { return _idx; }

            /// @returns edges of the node
            const auto &half_edges() const { return _half_edges; }

            /// Node(volume) index
            dindex _idx;
            /// Vector of half edges towards other volumes
            vector_t<typename detector_t::surface_type::mask_link> _half_edges;
        };

        /// @brief Iterator over the graph nodes.
        struct iterator {
            // Iterator type defs
            using value_type = node;
            using volume_iter =
                decltype(std::begin(std::declval<volume_container_t>()));

            /// Constructor from an iterator on the detector volume container
            /// and a reference to its surface container.
            iterator(volume_iter &&vol_itr, const surface_container_t &surfaces)
                : _vol_itr(vol_itr), _surfaces(surfaces) {}

            /// Dereference operator @returns a graph node
            node operator*() { return node(*_vol_itr, _surfaces); }

            /// Prefix increment. No postfix increment implemented
            iterator &operator++() {
                ++_vol_itr;
                return *this;
            }

            /// Equality operator
            bool operator==(const iterator &rhs) {
                return _vol_itr == rhs._vol_itr and
                       &_surfaces == &rhs._surfaces;
            }

            /// Inequality operator
            bool operator!=(const iterator &rhs) {
                return _vol_itr != rhs._vol_itr or &_surfaces != &rhs._surfaces;
            }

            /// Iterator over the detector volume container
            volume_iter _vol_itr;
            /// Detector surfaces
            const surface_container_t &_surfaces;
        };

        /// Collection value type
        using value_type = node;

        /// Constructor from a detector containers
        node_collection(const volume_container_t &volumes,
                        const surface_container_t &surfaces)
            : _volumes(volumes), _surfaces(surfaces) {}

        /// @returns the number of nodes
        std::size_t size() const { return _volumes.size(); }

        /// @returns the first node in the collection
        node front() const { return node(_volumes.front(), _surfaces); }

        /// @returns beginning of collection
        iterator begin() const { return iterator(_volumes.begin(), _surfaces); }

        /// @returns end of collection
        iterator end() const { return iterator(_volumes.end(), _surfaces); }

        const volume_container_t &_volumes;
        const surface_container_t &_surfaces;
    };

    /// @brief Builds graph edges from the detector mask collection on the fly.
    /// /// Iterates through the detectors mask store given a volume id (graph
    /// node) The graph edges are constructed for the node on the fly from the
    /// half edges in the @c _edges collection (detector masks).
    struct edge_collection {
        /// The link to the next volume
        using mask_edge_t = typename detector_t::surface_type::edge_type;
        /// The link to the surfaces mask
        using mask_link_t = typename detector_t::surface_type::mask_link;

        /// An edge in the graph connects two nodes (volumes). It is constructed
        /// from a surface and its mask.
        struct edge {
            edge(const dindex volume_id, const dindex link)
                : _from(volume_id), _to(link) {}

            dindex from() const { return _from; }
            dindex to() const { return _to; }

            dindex _from, _to;
        };

        /// @brief Iterator that constructs graph edges from mask links on the
        /// fly
        struct iterator {
            // Iterator type defs
            using value_type = edge;
            using edge_iter =
                decltype(std::begin(std::declval<vector_t<edge>>()));

            /// Constructor from a volume id and the masks of one of its
            /// surfaces
            iterator(const dindex volume_id, const mask_link_t &mask_link,
                     const edge_collection &edges) {
                build_edges_vector(volume_id, mask_link, edges.get_container());
                _itr = _edges.begin();
            }

            /// Dereference operator @returns a graph edge
            edge operator*() { return *_itr; }

            /// Prefix increment. No postfix increment implemented
            iterator &operator++() {
                ++_itr;
                return *this;
            }

            /// @returns begging of graph edges container
            iterator &begin() {
                _itr = _edges.begin();
                return *this;
            }

            /// @returns end of graph edges container
            iterator &end() {
                _itr = _edges.end();
                return *this;
            }

            /// Equality operator
            bool operator==(const iterator &rhs) {
                return _itr == rhs._itr and &_edges == &rhs._edges;
            }

            /// Inquality operator
            bool operator!=(const iterator &rhs) {
                return _itr != rhs._itr and &_edges != &rhs._edges;
            }

            /// @brief BUilds the collection of graph edges for a given node.
            ///
            /// From the volume index, a mask link owned by one of the volumes
            /// surfaces and the detector mask container, a vector of graph
            /// edges is built. The mask container needs to be unrolled to find
            /// the correct instance and therefore link to the next volume.
            ///
            /// @param volume_id the index of the volume/node
            /// @param mask_link a mask link of a surface belonging to that vol.
            /// @param masks the mask store of the detector
            template <std::size_t current_id = 0>
            inline void build_edges_vector(const dindex volume_id,
                                           const mask_link_t &mask_link,
                                           const mask_container_t &masks) {

                if (detail::get<0>(mask_link) == current_id) {
                    // Get the mask group
                    const auto &mask_group =
                        masks.template group<mask_container_t::to_id(
                            current_id)>();
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

            /// Iterator over the edges vector
            edge_iter _itr;
            /// Vector of graph edges, constructed on the fly in the iterator
            vector_t<edge> _edges;
        };

        /// Collection value type
        using value_type = edge;

        /// Constructor from the detector masks store.
        edge_collection(const typename detector_t::mask_container &masks)
            : _edges(masks) {}

        /// @return the container of half edges
        const mask_container_t &get_container() const { return _edges; }

        const mask_container_t &_edges;
    };

    // Graph nodes
    using node_type = typename node_collection::value_type;
    using node_iter = typename node_collection::iterator;
    // Graph edges
    using edge_type = typename edge_collection::value_type;
    using edge_iter = typename edge_collection::iterator;

    /// Default constructor
    volume_graph() = delete;

    /// @brief Build from existing nodes and edges, which are provide by the
    /// detector.
    ///
    /// @param det provides: geometry volumes that become the graph nodes,
    /// surfaces which are needed to index the correct masks and the
    /// masks that link to volumes and become graph edges.
    volume_graph(const detector_t &det)
        : _nodes(det.volumes(), det.surfaces()),
          _edges(det.mask_store()),
          _adj_matrix{0} {
        build();
    }

    /// Default destructor: we don't own anything.
    ~volume_graph() = default;

    /// @return number of nodes
    size_t n_nodes() const { return _nodes.size(); }

    /// @return node collection - const access. */
    const auto &nodes() const { return _nodes; }

    /// @return number of surfaces/portals in the geometry */
    // size_t n_edges() const { return _edges.size(); }

    /// @return edges collection - const access.
    const auto &edges() const { return _edges; }

    /// @return graph adjacency - const access.
    const auto &adjacency_matrix() const { return _adj_matrix; }

    /// Walks breadth first through the geometry objects.
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

    /// @returns the linking description as a string.
    inline const std::string to_string() const {
        std::stringstream stream;
        std::size_t dim = n_nodes() + 1;
        for (const auto &n : _nodes) {
            stream << "[>>] Node with index " << n.index() << std::endl;
            stream << " -> edges: " << std::endl;
            for (std::size_t i = 0; i < dim; ++i) {
                const auto degr = _adj_matrix[dim * n.index() + i];
                if (degr == 0) {
                    continue;
                }
                std::string n_occur =
                    degr > 1 ? "\t\t\t\t(" + std::to_string(degr) + "x)" : "";

                // Edge that leads out of the detector world
                if (i == dim - 1 and degr != 0) {
                    stream << "    -> leaving world " + n_occur << std::endl;
                } else {
                    stream << "    -> " << std::to_string(i) + "\t" + n_occur
                           << std::endl;
                }
            }
        }
        return stream.str();
    }

    private:
    /// @brief Go through the nodes and fill adjacency matrix.
    ///
    /// Root node is always at zero.
    void build() {
        // Leave space for the world volume (links to dindex_invalid)
        const std::size_t dim = n_nodes() + 1;
        _adj_matrix.resize(dim * dim);

        for (const auto &n : _nodes) {
            // Only works for non batched geometries
            for (const auto &edg_link : n.half_edges()) {
                // Build an edge
                for (const auto edg : edge_iter(n.index(), edg_link, _edges)) {
                    const dindex elem = edg.to() < dindex_invalid
                                            ? dim * n.index() + edg.to()
                                            : dim * n.index() + dim - 1;
                    _adj_matrix[elem]++;
                }
            }
        }
    }

    /// Graph nodes
    node_collection _nodes;

    /// Graph edges
    edge_collection _edges;

    /// Adjacency matrix
    vector_t<dindex> _adj_matrix = {};
};

}  // namespace detray