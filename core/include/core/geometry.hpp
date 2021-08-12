/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/surface_base.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>

namespace detray
{
    /**
     * @brief Indexed graph structure geometry.
     *
     * This class provides a geometry that treats volumes as graph nodes and
     * portals as edges. It contains a vector_type of volumes and portals, as
     * well as a vector of surfaces that either belong to a portal or are
     * surfaces contained in a volume.
     * All data in this geometry are either ranges or indices into the data
     * structures of the detector.
     *
     * @tparam array_type the type of the internal array, must have STL
     *                    semantics
     * @tparam vector_type the type of the internal array, must have STL
     *                     semantics
     * @tparam surface_source_link the type of the link to an external surface
     *                             source
     *
     * @note The geometry knows nothing about coordinate systems. This is
     *       handeled by dedicated geometry access objects.
     */
    template <template <typename, unsigned int> class array_type = darray,
              template <typename> class vector_type = dvector,
              typename surface_source_link = dindex>
    class index_graph_geometry
    {

    public:

        // In case the geometry needs to be printed
        using name_map = std::map<dindex, std::string>;
        using source_link = dindex;

        // Index container used in volumes {first index, number of elements}
        using element_range = array_type<dindex, 2>;

        /**
         * Half-edge in a volume graph with surface finders.
         *
         * \var next_volume points to adjacent volume by its index
         * \var next_obj_finder points to the correct object_finder in the
         *                      adjacent volume
         */
        struct half_edge
        {
            dindex next_volume     = dindex_invalid;
            dindex next_obj_finder = dindex_invalid;
        };

        /**
         * @brief The Surface definition
         *
         * Contains indices in the detector scope containers that define a
         * batch of surfaces.
         */
        struct surface_batch
        {
            /** Number of surfaces contained in this batch */
            dindex n_surfaces = dindex_invalid;
            /** Type of the masks for these surfaces */
            dindex mask_type  = dindex_invalid;
            /**
              * Start of the surface masks in the detector scope container
              * and number of masks per surface.
              */
            element_range mask_range = {dindex_invalid, 0};
            /**
              * Start of the surface transforms in the detector scope
              * container
              */
            dindex transform_idx = dindex_invalid;
            /** Start of the surface links in the geometry owned container */
            dindex source_idx    = dindex_invalid;

            /**
             * Get the mask range in the detector container for a single surface
             * in the batch.
             *
             * @param surface_idx the surface idx relative to the batch
             *
             * @return the mask range in the detector container.
             */
            element_range mask_range_by_surface(const dindex &surface_index)
            const
            {
                dindex first_mask = mask_range[0]
                                    + mask_range[1] * surface_index;
                return {first_mask, first_mask + mask_range[1]};
            }
        };

        /**
         * The Portal is an edge with an associated surface between two
         * volumes. This struct represents one half of that edge, grouped by
         * mask type.
         *
         * \var half-edge_batch ponits to a range of half-edges for this mask
         *                      type
         * \var surface_batch ponits to a corresponding surface batch in the
         *                    surfaces container
         */
        struct portal_batch
        {
            /**
             * Start index in the geometry owned edges container. The edges
             * are sorted according to the surface masks.
             */
            dindex half_edges = dindex_invalid;
            /** 
             * The surface batches asociated with this portal. One batch per
             * mask type. The number of portal links is included in the
             * surface batch.
             */
            element_range surface_batches = {dindex_invalid, 0};
        };

        /**
         * Nested volume class that holds the local information of a
         * volume, its portals and the contained surfaces.
         */
        struct volume
        {
            /** Deleted constructor */
            volume() = delete;

            // Allowed constructors

            /**
             * Contructor with name and bounds
             *
             * @param bounds of the volume
             */
            volume(const array_type<scalar, 6> &bounds) : _bounds(bounds) {}

            /** @return the bounds - const access */
            const array_type<scalar, 6> &bounds() const { return _bounds; }

            /** @return the volume index */
            dindex index() const { return _index; }

            /** @return the entry into the local surface finders */
            const dindex &surfaces_finder_entry() const { return _surfaces_finder; }

            /** @return whether the volume is empty or not */
            bool empty() const { return _surfaces[1] == dindex_invalid; }

            /** @return number of surfaces */
            const auto &n_surface_batches() const { return _surfaces[1]; }

            /** @return range in surfaces contianer */
            const auto &surface_range() const { return _surfaces; }

            /** @return number of portals */
            const auto &n_portal_batches() const { return _portals.surface_batches[1]; }

            /** @return range in surfaces contianer for volume portals */
            const auto &portal_range() const { return _portals.surface_batches; }

            /** @return the portals themselves */
            const auto &portals() const { return _portals; }

            /**
             * Set the index into the surfaces container.
             *
             * @param surface_range Range of surfaces in surface container
             */
            void add_surfaces(element_range surface_range)
            {
                // initialize volume
                if (_surfaces[0] == dindex_invalid)
                {
                    _surfaces[0] = surface_range[0];
                    _surfaces[1] = surface_range[1];
                }
                // volume only keeps track of its range
                else
                {
                    _surfaces[1] += surface_range[1];
                }
            }

            /**
             * Set the index into the geometry and detector data structures
             * for portals
             *
             * @param portal_range Range of portal structs in portal container
             */
            void add_portals(portal_batch portals)
            {
                // initialize volume
                if (_portals.half_edges == dindex_invalid)
                {
                    _portals = portals;
                }
                // portal batch range has widened, but all initial indices are 
                // the same
                else
                {
                    _portals.surface_batches[1] += portals.surface_batches[1];
                }
            }

            /**
             * Set the index into the detector transform store for portals
             *
             * @param finder_index Index of the local object finder
             */
            void add_surfaces_finder(dindex finder_index)
            {
                _surfaces_finder = finder_index;
            }

            /** Bounds section, default for r, z, phi */
            array_type<scalar, 6> _bounds =
            {
                0.,
                std::numeric_limits<scalar>::max(),
                -std::numeric_limits<scalar>::max(),
                std::numeric_limits<scalar>::max(),
                -M_PI, M_PI
            };

            /** Volume index: Position in geometry volume container */
            dindex _index = dindex_invalid;

            /** Range of indices for surfaces in surfaces container */
            element_range _surfaces = {dindex_invalid, 0};

            /** Range of indices for portals in portals container */
            portal_batch _portals = {};

            /** Index into the surface finder container */
            dindex _surfaces_finder = dindex_invalid;
        };

        /**
         * Add a new volume and retrieve a reference to it.
         *
         * @param bounds of the volume, they are expected to be already
         *               attaching
         * @param surfaces_finder_entry of the volume, where to entry the
         *                              surface finder
         *
         * @return non-const reference of the new volume
         */
        volume &new_volume(const array_type<scalar, 6> &bounds,
            dindex surfaces_finder_entry = dindex_invalid)
        {
            _nodes.push_back(std::move(volume(bounds)));
            dindex cvolume_idx = _nodes.size() - 1;
            volume &cvolume    = _nodes[cvolume_idx];
            cvolume._index     = cvolume_idx;
            cvolume._surfaces_finder = surfaces_finder_entry;

            return cvolume;
        }

        /** @return total number of nodes (volumes) */
        const size_t n_volumes()  const { return _nodes.size(); }

        /** @return all volumes in the geometry - const access. */
        const vector_type<volume> &volumes() const { return _nodes; }

        /** @return the volume by @param volume_index - const access. */
        const volume &volume_by_index(dindex volume_index) const { return _nodes[volume_index]; }

        /** @return the volume by @param volume_index - non-const access. */
        volume &volume_by_index(dindex volume_index) { return _nodes[volume_index]; }

        /** @return number of surface batches - const access. */
        const auto n_surface_batches() const { return _surfaces.size(); }

        /** @return number of surfaces in batch range. */
        const dindex n_surfaces_in_range(const element_range &range) const
        {
            // Get the number of surfaces in every batch
            dindex n_sfs = 0;
            for (size_t bi = range[0]; bi < range[0] + range[1]; bi++)
            {
                n_sfs += _surfaces[bi].n_surfaces;
            }
            return n_sfs;
        }

        /** @return number of surfaces */
        const dindex n_surfaces(const volume &v) const { return n_surfaces_in_range(v.surface_range()); }

        /** @return number of surfaces */
        const dindex n_portals(const volume &v) const { return n_surfaces_in_range(v.portal_range()); }

        /** @return all surface batches in geometry. */
        const auto &surfaces() const { return _surfaces; }

        /** @return all graph edges - const access */
        const auto &edges() const { return _edges; }

        /** @return all source links - const access */
        const auto &source_links() const { return _source_links; }

        /**
         * Find the geometry graph edge for a given a certain surface and
         * mask.
         *
         * @param index The index of the surface in the portal batch
         *
         * @return the graph (half) edge belonging to the index.
         */
        const inline auto get_edge(const portal_batch &pb,
                                   const std::array<dindex, 3> &index) const
        {
            // Get the number of masks in the previous surface batches
            // in this portal
            dindex mask_index = 0;
            for (size_t bi = pb.surface_batches[0]; bi < index[0]; bi++)
            {
                mask_index += _surfaces[bi].n_surfaces *
                              _surfaces[bi].mask_range[1];
            }
            // Add the mask offset corresponding to the passed index
            const auto &batch = _surfaces[index[0]];
            mask_index += batch.mask_range[1] * index[1] + index[2];
            // Return the graph edge corresponding to this mask
            return _edges[pb.half_edges + mask_index];
        }

        /**
         * Find the next volume for a given a certain surface and mask.
         *
         * @param index The index of the surface in the portal batch
         *
         * @return the next volume index according to the index.
         */
        const inline dindex next_volume(const portal_batch pb,
                                        const std::array<dindex, 3> index) 
            const { return get_edge(pb, index).next_volume; }

        /**
         * Find the next surfaces finder for a given a certain surface and
         * mask.
         *
         * @param index The index of the surface in the portal batch
         *
         * @return the next volume index according to the index.
         */
        const inline dindex next_surfaces_finder(const portal_batch pb,
                                            const std::array<dindex, 3> index)
        const { return get_edge(pb, index).next_obj_finder; }

        /**
          * Adds a number of detector/portal surfaces to a specific volume.
          *
          * @param volume The volume to which the new surfaces belong
          * @param surfaces The new surfaces
          * @param source_links All source links belonging to the surfaces
          */
        void add_surfaces(volume& volume, vector_type<surface_batch>& surfaces,
                          vector_type<surface_source_link>& source_links,
                          bool is_portal_surfaces = false)
        {
            auto sf_start = _surfaces.size();
            auto sl_start = _source_links.size();

            _surfaces.reserve(sf_start + surfaces.size());
            _surfaces.insert(_surfaces.end(), surfaces.begin(), surfaces.end());

            _source_links.reserve(sl_start + source_links.size());
            _source_links.insert(_source_links.end(), source_links.begin(), source_links.end());

            // Update source link start index to global container
            for (auto &sf : surfaces)
            {
                sf.source_idx += sl_start;
            }

            // Portal surfaces are recorded in the portal batch
            if (not is_portal_surfaces)
            {
                volume.add_surfaces({sf_start, surfaces.size()});
            }
        }

        /**
          * Adds a number of portals to a specific volume. The surface links
          * of portals need to be global surface links in the surfaces
          * container. The mask links have been set correctly by the detector.
          * 
          * @param volume The volume to which the new surfaces belong
          * @param portals The new portals: A half-edge per surface in the
          *                batch
          * @param edges The half-edges for every surface in the batch
          */
        void add_portals(volume& volume, portal_batch& portals,
                         vector_type<half_edge>& edges)
        {
            auto edg_start = _edges.size();

            _edges.reserve(edg_start + edges.size());
            _edges.insert(_edges.end(),  edges.begin(), edges.end());

            portals.half_edges += edg_start;
            
            volume.add_portals(portals);
        }

        /**
          * Add a surface finder index to a volume.
          *
          * @param volume The volume to which the index should be added.
          * @param finder_index The index of the surface finder for the volume.
          */
        void add_surfaces_finder(volume& volume, dindex& finder_index)
        {
            volume.add_surfaces_finder(finder_index);
        }

        /**
          * Print geometry if an external name map is provided for the volumes.
          *
          * @param names  Lookup for the names by volume index.
          *
          * @returns the geometry description as a string
          */
        template <typename name_map>
        const std::string to_string(name_map &names) const
        {
            std::stringstream ss;
            for (const auto &[i, v] : enumerate(_nodes))
            {
                ss << "[>>] Volume at index " << i
                   << " - name: '"   << names[v.index()]
                   << "'" << std::endl;

                ss << "     contains    "   << v.n_surface_batches()
                   << " surface batches ("  << n_surfaces(v)
                   << " detector surfaces)" << std::endl;

                ss << "                 "         << v.n_portal_batches()
                   << " portal surface batches (" << n_portals(v)
                   << " detector portals)"        << std::endl;
                   
                if (v.surfaces_finder_entry() != dindex_invalid)
                {
                    ss << "  sf finders idx " << v.surfaces_finder_entry()
                       << std::endl;
                }
                ss << "     bounds r = (" << v._bounds[0] << ", "
                   << v._bounds[1] << ")" << std::endl;
                ss << "            z = (" << v._bounds[2] << ", "
                   << v._bounds[3] << ")" << std::endl;
            }
            return ss.str();
        };

    private:

        /** Surface objects, including those for portal surfaces. */
        vector_type<surface_batch> _surfaces = {};

        /** Source_links for surfaces. */
        vector_type<surface_source_link> _source_links = {};

        /** The half-edges of the geometry graph (volume portals). */
        vector_type<half_edge> _edges = {};

        /** Contains the geometrical relations encoded in volume nodes. */
        vector_type<volume> _nodes = {};
    };

}
