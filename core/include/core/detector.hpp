/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/geometry.hpp"
#include "core/intersection.hpp"
#include "core/surface_base.hpp"
#include "core/transform_store.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "masks/masks.hpp"
#include "utils/indexing.hpp"
#include "utils/enumerate.hpp"
#include "tools/planar_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/local_object_finder.hpp"

#include <cassert>
#include <iterator>
#include <string>
#include <sstream>
#include <type_traits>

namespace detray
{

    // Algebra, point2 is not strongly typed
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using point2 = __plugin::point2;

    /** Indexed detector definition.
     * 
     * This class is a heavy templated detector definition class, that connects
     * surfaces, layers and volumes via an indexed system.
     *
     * @tparam array_type the type of the internal array, must have STL semantics
     * @tparam tuple_type the type of the internal tuple, must have STL semantics
     * @tparam vector_type the type of the internal array, must have STL semantics
     * @tparam alignable_store the type of the transform store
     * @tparam surface_source_link the type of the link to an external surface source
     * @tparam bounds_source_link the type of the link to an external bounds source
     * @tparam surfaces_populator_type the type of populator used to fill the surfaces grids
     * @tparam surfaces_serializer_type the type of the memory serializer for the surfaces grids
     * 
     */
    template <template <typename, unsigned int> class array_type = darray,
              template <typename...> class tuple_type = dtuple,
              template <typename> class vector_type = dvector,
              typename alignable_store = static_transform_store<vector_type>,
              typename surface_source_link = dindex,
              typename bounds_source_link = dindex,
              typename geometry_t = index_graph_geometry<array_type, vector_type, dindex>,
              typename surfaces_populator_type = attach_populator<false, dindex, vector_type>,
              typename surfaces_serializer_type = serializer2>
    class detector
    {

    public:
        /** Encodes the position in a collection container for the respective
            type. */
        enum mask_container_index : unsigned int {
            e_mask_types = 5,
            e_rectangle2 = 0,
            e_trapezoid2 = 1,
            e_annulus2 = 2,
            e_cylinder3 = 3,
            e_ring2 = 4,
            e_single3 = std::numeric_limits<unsigned int>::max(),
            e_unknown = std::numeric_limits<unsigned int>::max(),
        };

        /// Forward the alignable container and context
        using transform_store = alignable_store;
        using context = typename alignable_store::context;


        // Geometry defined types
        using portal  = typename geometry_t::portal_batch;
        using surface = typename geometry_t::surface_batch;
        using volume  = typename geometry_t::volume;
        using mask_index = typename geometry_t::element_range;
        // Link to an external source
        using surface_links = typename geometry_t::source_link;
        using portal_links = typename geometry_t::half_edge;


        /// Volume grid definition
        using volume_grid = grid2<replace_populator<dindex, std::numeric_limits<dindex>::max(), vector_type>,
                                  axis::irregular<array_type, vector_type>,
                                  axis::irregular<array_type, vector_type>,
                                  serializer2>;

        using surfaces_regular_axis = axis::regular<array_type>;
        using surfaces_circular_axis = axis::circular<array_type>;
        using surfaces_regular_circular_grid = grid2<surfaces_populator_type,
                                                     surfaces_regular_axis,
                                                     surfaces_circular_axis,
                                                     surfaces_serializer_type,
                                                     array_type,
                                                     tuple_type,
                                                     vector_type>;

        using surfaces_finder = surfaces_regular_circular_grid;


        // (Portal-)surface types (masks)
        using rectangle = rectangle2<planar_intersector, __plugin::cartesian2>;
        using trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2>;
        using annulus = annulus2<planar_intersector, __plugin::cartesian2>;
        using cylinder = cylinder3<false, cylinder_intersector, __plugin::cylindrical2>;
        using disc = ring2<planar_intersector, __plugin::cartesian2>;

        // Container definitions
        using link_container = vector_type<surface_links>;
        using mask_container = tuple_type<vector_type<rectangle>,
                                          vector_type<trapezoid>,
                                          vector_type<annulus>,
                                          vector_type<cylinder>,
                                          vector_type<disc>>;

        // Temporary containers (sorted by mask type) to fill the detector
        using portal_container = vector_type<portal_links>;
        using surface_container = array_type<surface, e_mask_types>;
        using transform_container = tuple_type<
                                        typename alignable_store::storage,
                                        typename alignable_store::storage,
                                        typename alignable_store::storage,
                                        typename alignable_store::storage,
                                        typename alignable_store::storage>;

        /** Allowed costructor
         * @param name the detector
         */
        detector(const std::string &name) : _name(name) {}

        /** Allowed costructor
         * @param name the detector
         */
        detector(const std::string &name, geometry_t geo)
            : _name(name), _geometry(std::move(geo)) {}

        /** Copy constructor makes sure the volumes belong to new detector.
         *
         * @param other Detector to be copied
         */
        detector(const detector &other) = default;
        detector() = delete;
        ~detector() = default;

        /** @return the name of the detector */
        const geometry_t &geometry() const { return _geometry; }

        /** @return the name of the detector */
        const std::string &name() const { return _name; }

        /** @return the contained volumes of the detector - const access */
        const vector_type<volume> &volumes() const { return _geometry.volumes(); }

        /** @return the volume by @param volume_index - const access */
        const volume &volume_by_index(dindex volume_index) const { return _geometry.volume_by_index(volume_index); }

        /** @return the volume by @param volume_index - non-const access */
        volume &volume_by_index(dindex volume_index) { return _geometry.volume_by_index(volume_index); }

        /** @return the volume by @param position - const access */
        const volume &volume_by_pos(const point3 &p) const
        {
            point2 p2 = {getter::perp(p), p[2]};
            dindex volume_index = _volume_grid.bin(p2);
            return _geometry.volume_by_index(volume_index);
        }

        /** @return all surfaces - const access */
        const auto &surfaces() const { return _geometry.surfaces(); }

        /** @return all surfaces - const access */
        const auto &portals() const { return _geometry.edges(); }

        /** @return all masks - const access */
        const auto &masks() const { return _masks; }

        /** @return all surfaces - const access */
        const auto &transforms(const context &/*ctx*/) const { return _transforms; }

        /** Get the current transform index for surfaces
         *
         * @param ctx The context of the call
         *
         * @return Index to add new transforms at 
         */
        const unsigned int transform_index(const context &ctx) const { return _transforms.size(ctx); }

        /** @return the volume grid - const access */
        const volume_grid &volume_search_grid() const { return _volume_grid; }

        /** @return the surface finders - const access */
        const vector_type<surfaces_finder> &surfaces_finders() const { return _surfaces_finders; }

        /** Add a new volume and retrieve a reference to it
         *
         * @param name of the volume
         * @param bounds of the volume, they are expected to be already
         *               attaching
         * @param surfaces_finder_entry of the volume, where to entry the
         *                              surface finder 
         *
         * @return non-const reference of the new volume
         */
        volume &new_volume(const array_type<scalar, 6> &bounds, dindex surfaces_finder_entry = dindex_invalid)
        {
            return _geometry.new_volume(bounds, surfaces_finder_entry);
        }

        /**
         * Add a new full set of alignable transforms for surfaces
         *
         * @param volume The volume the portals belong to
         * @param surfaces New surfaces (one batch per mask type)
         * @param masks Corresponding surface masks
         * @param trfs Corresponding surface transforms (same order as surfaces
         *             in batch!)
         * @param sc_links Corresponding surface source links (same order as
         *                 surfaces in batch!)
         * @param ctx Transform context
         *
         * @note can throw an exception if input data is inconsistent
         */
        void add_surfaces(
            volume &volume,
            surface_container &surfaces,
            mask_container &masks,
            transform_container &trfs,
            link_container &sc_links,
            const typename alignable_store::context ctx = {},
            bool is_portal_surfaces = false) noexcept(false)
        {
            // Unroll the intersection depending on the mask container size
            vector_type<surface> valid_surfaces;
            unroll_container_filling(surfaces, masks, trfs, valid_surfaces, ctx);

            _geometry.add_surfaces(volume, valid_surfaces, sc_links, is_portal_surfaces);
        }

        /**
         * Add a new full set of alignable transforms for surfaces.
         *
         * @param volume The volume the portals belong to
         * @param edges The connections between volumes
         * @param surfaces New portal surfaces that might need to be added
         * @param masks If we need new surfaces, add the masks
         * @param trfs If we need new surfaces, add the transforms
         * @param sc_links Links to external source for surfaces
         *
         * @note can throw an exception if input data is inconsistent
         */
        void add_portals(
            volume &volume,
            vector_type<typename geometry_t::half_edge> &edges,
            surface_container &surfaces,
            mask_container &masks,
            transform_container &trfs,
            link_container &sc_links) noexcept(false)
        {
            portal portals = { 0, {_geometry.n_surface_batches(),
                                   _geometry.n_surface_batches()} };

            // TODO: Put into geometry consistency check
            // Ensure that there is an edge for every mask
            unsigned int n_sfs = 0;
            for (const auto &sfb : surfaces) { n_sfs += sfb.mask_range[1]; }
            assert(n_sfs == edges.size());

            // TODO: Implement portal surface sharing
            if (not sc_links.empty())
            {
                add_surfaces(volume, surfaces, masks, trfs, sc_links, {}, true);
                // How many batches were added for the portals?
                portals.surface_batches[1] = _geometry.n_surface_batches() - portals.surface_batches[0];
            }

            _geometry.add_portals(volume, portals, edges);
        }

        /**
         * Add a new set of surface transforms to the store.
         *
         * @param ctx The current context for the transforms
         * @param trfs The transforms that need to be added
         *
         * @note can throw an exception if input data is inconsistent
         */
        void add_transforms(const context &ctx,
                            typename alignable_store::storage &trfs) noexcept(false)
        {
            _transforms.append_contextual_transforms(ctx, trfs);
        }

        /** Add the volume grid - move semantics 
         *
         * @param v_grid the volume grid to be added
         */
        void add_volume_grid(volume_grid &&v_grid) noexcept
        {
            _volume_grid = std::move(v_grid);
        }

        /**
         * Add local surface finders linked to from the portals
         * - move semantics
         *
         * This connects portals and surface grids
         */
        void add_surfaces_finders(vector_type<surfaces_finder> &&surfaces_finders) noexcept
        {
            _surfaces_finders = std::move(surfaces_finders);
        }

        /** @return Output to string */
        template<typename name_map>
        const std::string to_string(name_map &names) const
        {
            std::stringstream ss;
            ss << "[>] Detector '" << _name << "' has "
               << _geometry.n_volumes() << " volumes," << std::endl;

            ss << "    contains  " << _surfaces_finders.size()
               << " local surface finders." << std::endl;

            ss << _geometry.to_string(names) << std::endl;
            return ss.str();
        };

    private:

        template <size_t current_idx = 0>
        void unroll_container_filling(surface_container &surfaces,
                                      mask_container &masks,
                                      transform_container &trfs,
                                      vector_type<surface> &valid_surfaces,
                                      const typename alignable_store::context ctx = {})
        {
            auto surface_batch = surfaces[current_idx];
            // Skip uninitialized batches and empty batches
            if (surface_batch.n_surfaces != dindex_invalid)
            {
                const auto &surface_masks      = std::get<current_idx>(masks);
                auto &detector_masks           = std::get<current_idx>(_masks);
                const auto &surface_transforms = std::get<current_idx>(trfs);

                if (surface_transforms.size() != 0 and surface_masks.size() != 0)
                {
                    // Fill into detectors containers
                    const auto first_mask = detector_masks.size();
                    const auto first_trsf = transform_index(ctx);

                    detector_masks.reserve(first_mask + surface_masks.size());
                    detector_masks.insert(detector_masks.end(), surface_masks.begin(),
                                        surface_masks.end());

                    _transforms.append_contextual_transforms(ctx, surface_transforms);

                    // Update batch
                    surface_batch.n_surfaces = surface_transforms.size();
                    surface_batch.mask_type = current_idx;
                    surface_batch.mask_range[0] = surface_transforms.size() == 0 ?
                                                dindex_invalid : surface_batch.mask_range[0] + first_mask;
                    surface_batch.transform_idx = surface_transforms.size() == 0 ?
                                                dindex_invalid : surface_batch.transform_idx + first_trsf;

                    valid_surfaces.push_back(surface_batch);
                }
            }
            // Next mask type
            if constexpr (current_idx < std::tuple_size_v<mask_container> - 1) {
                return unroll_container_filling<current_idx + 1> (surfaces, masks, trfs, valid_surfaces, ctx);
            }
        }

        std::string _name = "unknown";

        /** Geometry struct keeps the relations between geometry objects */
        geometry_t _geometry{};

        /** Keeps all of the surfaces transform data in contiguous memory */
        transform_store _transforms{};

        /** Contains the masks for all surfaces, sorted by type */
        mask_container _masks;

        vector_type<surfaces_finder> _surfaces_finders;

        volume_grid _volume_grid = volume_grid(std::move(axis::irregular{{}}),
                                               std::move(axis::irregular{{}}));
    };

} // namespace detray
