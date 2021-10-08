/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <sstream>
#include <string>
#include <utility>

#include "geometry/unified_index_geometry.hpp"
#include "masks/masks.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {
/**
 * @brief Uses the geometry implementations to walk through their graph-like
 * structure breadth first.
 *
 * This class provides a geometry that defines logic volumes which contain
 * the detector surfaces, joined together by dedicated portal surfaces. It
 * exports all types needed for navigation and strictly only keeps the
 * index data (links) that define the geometry relations. The simple geometry
 * itself makes no distinction between surfaces and portals. Both carry the
 * same link type: a portal points to the next volume, a surface to the current
 * volume
 *
 * @tparam array_type the type of the internal array, must have STL
 *                    semantics
 * @tparam vector_type the type of the internal array, must have STL
 *                     semantics
 * @tparam surface_source_link the type of the link to an external surface
 *                             source
 *
 * @note The geometry knows nothing about coordinate systems. This is
 *       handeled by geometry access objects (e.g. the grid).
 */
/*template <template <typename, unsigned int> class array_type = darray,
          template <typename> class vector_type = dvector,
          template <typename...> class tuple_type = dtuple,
          typename geometry_type = simple_geometry<array_type, vector_type,
                                                  tuple_type, dindex, dindex>>*/
template<typename geometry,
         template <typename> class vector_type = dvector>
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
    geometry_graph(const vector_type<typename geometry::volume_type> &volumes, const vector_type<typename geometry::portal> &portals): _nodes(volumes), _edges(portals){}

    /** Default destructor: we don't own anything */
    ~geometry_graph() = default;

    /** @return total number of volumes */
    const size_t n_nodes() const { return _nodes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto nodes() const { return _nodes; }

    /** @return all surfaces/portals in the geometry */
    size_t n_edges() const {
        return _edges.size();
    }

    /** @return all surfaces/portals in the geometry */
    const auto edges() const {
        return _edges;
    }

    /**
     * Print geometry if an external name map is provided for the volumes.
     *
     * @param names  Lookup for the names by volume index.
     *
     * @returns the geometry description as a string
     */
    inline const std::string to_string() const {
        std::stringstream ss;
        for (const auto &n : enumerate(_nodes)) {
            ss << "[>>] Node with index " << n.index() << std::endl;
        }
        return ss.str();
    };

    private:

    /** Graph nodes */
    const vector_type<node_t>& _nodes;

    /** Graph edges */
    const vector_type<edge_t>& _edges;
};

}  // namespace detray