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
#include "detray/utils/ranges.hpp"

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
    struct node_generator
        : public detray::ranges::view_interface<node_generator> {

        /// A node in the graph has an index (volume index) and a collction of
        /// edge that belong to it (mask link of every surface in the volume).
        /// One mask link is a half edge in the graph.
        struct node {

            /// Constructor from a detectors volume and surface collections
            node(const typename detector_t::volume_type &volume,
                 const surface_container_t &surfaces)
                : _idx(volume.index()) {
                for (const auto &sf :
                     detray::ranges::subrange(surfaces, volume)) {
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
            using volume_iter =
                detray::ranges::const_iterator_t<volume_container_t>;

            using difference_type =
                typename std::iterator_traits<volume_iter>::difference_type;
            using value_type = node;
            using pointer = typename std::iterator_traits<volume_iter>::pointer;
            using reference =
                typename std::iterator_traits<volume_iter>::reference;
            using iterator_category =
                typename std::iterator_traits<volume_iter>::iterator_category;

            /// Constructor from an iterator on the detector volume container
            /// and a reference to its surface container.
            iterator(volume_iter &&vol_itr, const surface_container_t &surfaces)
                : _vol_itr(vol_itr), _surfaces(surfaces) {}

            /// Equality operator
            bool operator==(const iterator &rhs) const {
                return _vol_itr == rhs._vol_itr and
                       &_surfaces == &rhs._surfaces;
            }

            /// Inequality operator
            bool operator!=(const iterator &rhs) const {
                return not(*this == rhs);
            }

            /// Dereference operator @returns a graph node
            node operator*() { return node(*_vol_itr, _surfaces); }

            /// Prefix increment. No postfix increment implemented
            iterator &operator++() {
                ++_vol_itr;
                return *this;
            }

            /// Prefix decrement. No postfix decrement implemented
            iterator &operator--() {
                --_vol_itr;
                return *this;
            }

            /// @returns an iterator that has been advanced by @param j
            constexpr auto operator+(const difference_type j) const
                -> iterator {
                return {_vol_itr + j, _surfaces};
            }

            /// @returns an iterator that has been advanced by - @param j
            constexpr auto operator-(const difference_type j) const
                -> iterator {
                return *this + -j;
            }

            /// @returns distance between two iterators
            constexpr auto operator-(const iterator &other) const
                -> difference_type {
                return _vol_itr - other._vol_itr;
            }

            /// Advances iterator by @param j
            constexpr auto operator+=(const difference_type j) -> iterator & {
                _vol_itr += j;
                return *this;
            }

            /// Advances iterator by - @param j
            constexpr auto operator-=(const difference_type j) -> iterator & {
                return *this += -j;
            }

            /// Iterator over the detector volume container
            volume_iter _vol_itr;
            /// Detector surfaces
            const surface_container_t &_surfaces;
        };

        /// Node iterator type
        using iterator_t = iterator;

        /// No default constructor (holds member reference)
        node_generator() = delete;

        /// Constructor from a detector containers
        node_generator(const volume_container_t &volumes,
                       const surface_container_t &surfaces)
            : _volumes(volumes), _surfaces(surfaces) {}

        /// @returns beginning of collection
        iterator begin() const { return iterator(_volumes.begin(), _surfaces); }

        /// @returns end of collection
        iterator end() const { return iterator(_volumes.end(), _surfaces); }

        /// @returns the number of nodes
        std::size_t size() const { return _volumes.size(); }

        const volume_container_t &_volumes;
        const surface_container_t &_surfaces;
    };

    /// @brief Builds graph edges from the detector mask collection on the fly.
    ///
    /// Iterates through the detectors mask store given a volume id (graph
    /// node) The graph edges are constructed for the node on the fly from the
    /// half edges in the @c _edges collection (detector masks).
    struct edge_generator
        : public detray::ranges::view_interface<edge_generator> {
        /// The link to the next volume
        using mask_edge_t = typename detector_t::surface_type::volume_link_type;
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

        /// Nested functor that fills the edges from a mask container
        struct edges_builder {

            using output_type = vector_t<edge>;

            /// @brief Builds the collection of graph edges for a given
            /// node.
            ///
            /// From the volume index and a mask link owned by one of the
            /// volumes surfaces a vector of graph edges is built.
            /// The mask container calls this functor and provides the correct
            /// mask group from which the surface masks can be obtained from
            /// the surfaces mask range.
            ///
            /// @param mask_group the group of masks in the mask container
            ///                   of the detector
            /// @param volume_id the index of the volume/node for which to
            ///                  build edges
            /// @param mask_range the range of masks in the group for which
            ///                   to build edges
            template <typename mask_group_t, typename mask_range_t>
            inline output_type operator()(const mask_group_t &mask_group,
                                          const mask_range_t &mask_range,
                                          const dindex volume_id) {
                vector_t<edge> edges{};
                for (const auto &mask :
                     detray::ranges::subrange(mask_group, mask_range)) {
                    edges.emplace_back(volume_id, mask.volume_link());
                }
                return edges;
            }
        };

        /// No default constructor (holds member reference)
        edge_generator() = delete;

        /// Constructor from the detector masks store.
        edge_generator(const typename detector_t::mask_container &masks,
                       const dindex volume_id = 0,
                       const mask_link_t mask_link = {})
            : _masks(masks), _volume_id(volume_id), _mask_link{mask_link} {
            _edges = _masks.template call<edges_builder>(
                _mask_link, _volume_id);
        }

        /// @returns begging of graph edges container
        auto begin() { return _edges.begin(); }

        /// @returns end of graph edges container
        auto end() { return _edges.end(); }

        /// @return updated version of @c edge_generator
        edge_generator &operator()(dindex volume_id,
                                   const mask_link_t mask_link) {
            _volume_id = volume_id;
            _mask_link = mask_link;
            _edges.clear();
            _edges = _masks.template call<edges_builder>(
                _mask_link, _volume_id);

            return *this;
        }

        const mask_container_t &_masks;
        vector_t<edge> _edges;
        dindex _volume_id;
        mask_link_t _mask_link;
    };

    // Graph nodes
    using node_type = typename node_generator::node;
    // Graph edges
    using edge_type = typename edge_generator::edge;

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
    size_t n_edges(dindex volume_id) const {
        const node_type n = _nodes[volume_id];
        // Number of half_edges is equal to number of edges for volume
        return n.half_edges().size();
    }

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
    edge_generator::iterator(current->index(), edg_link, _edges)) { dindex nbr
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
                for (const auto edg : _edges(n.index(), edg_link)) {
                    const dindex elem = edg.to() < dindex_invalid
                                            ? dim * n.index() + edg.to()
                                            : dim * n.index() + dim - 1;
                    _adj_matrix[elem]++;
                }
            }
        }
    }

    /// Graph nodes
    node_generator _nodes;

    /// Graph edges
    edge_generator _edges;

    /// Adjacency matrix
    vector_t<dindex> _adj_matrix = {};
};

}  // namespace detray