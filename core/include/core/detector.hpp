/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

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

#include <string>
#include <sstream>

namespace detray
{

    // Algebra, point2 is not strongly typed
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using point2 = __plugin::point2;

    /** Indexed detector definition.
     *
     * @tparam alignable_store the type of the transform store
     * @tparam surface_source_link the type of the link to an external surface source
     * @tparam bounds_source_link the type of the link to an external bounds source
     * 
     */
    template <typename alignable_store = static_transform_store<>,
              typename surface_source_link = dindex,
              typename bounds_source_link = dindex,
              template <typename, unsigned int> class array_type = darray,
              template <typename...> class tuple_type = dtuple,
              template <typename> class vector_type = dvector>
    class detector
    {

    public:
        /// Forward the alignable container and context
        using transform_store = alignable_store;
        using context = typename alignable_store::context;

        /// Volume grid definition
        using volume_grid = grid2<replace_populator<>, axis::irregular<>, axis::irregular<>, serializer2>;

        /// Portals components:
        /// - links:  next volume, next (local) object finder
        using portal_links = array_type<dindex, 2>;
        /// - masks, with mask identifiers 0, 1
        using portal_cylinder = cylinder3<false, cylinder_intersector, __plugin::cylindrical2, portal_links, 0>;
        using portal_disc = ring2<planar_intersector, __plugin::cartesian2, portal_links, 1>;
        // - mask index: type, { first/last }
        using portal_mask_index = tuple_type<dindex, array_type<dindex, 2>>;
        using portal_mask_container = tuple_type<vector_type<portal_cylinder>, vector_type<portal_disc>>;

        /** The Portal definition:
         *  <transform_link, mask_index, volume_link, source_link >
         * 
         * transform_link: index into the transform container
         * mask_index: typed index into the mask container
         * volume_link: index of the volume this portal belongs to
         * source_link: some link to an eventual exernal representation
         * 
         */
        using portal = surface_base<dindex, portal_mask_index, dindex, surface_source_link>;
        using portal_container = vector_type<portal>;

        /// Surface components:
        /// - surface links
        using surface_links = array_type<dindex, 1>;
        /// - masks, with mask identifiers 0,1,2
        using surface_rectangle = rectangle2<planar_intersector, __plugin::cartesian2, surface_links, 0>;
        using surface_trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, surface_links, 1>;
        using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2, surface_links, 2>;
        using surface_cylinder = cylinder3<false, cylinder_intersector, __plugin::cylindrical2, surface_links, 3>;
        /// - mask index: type, entry
        using surface_mask_index = array_type<dindex, 2>;
        using surface_mask_container = tuple_type<vector_type<surface_rectangle>,
                                                  vector_type<surface_trapezoid>,
                                                  vector_type<surface_annulus>,
                                                  vector_type<surface_cylinder>>;

        using surface_link = surface_source_link;
        /** The Surface definition:
         *  <transform_link, mask_link, volume_link, source_link >
         */
        using surface = surface_base<dindex, surface_mask_index, dindex, surface_link>;
        using surface_container = vector_type<surface>;

        using surface_neighborhood = array_type<array_type<scalar, 2>, 2>;
        using surface_finder = std::function<vector_type<dindex>(const point2 &, const surface_neighborhood &nhood)>;

        /** Nested volume struct that holds the local information of the
         * volume and its portals.
         */
        class volume
        {

            friend class detector<alignable_store, surface_source_link, bounds_source_link>;

        public:
            /** Object holder class to synchronize access 
             * 
             * @tparam object_t the object type, either surface or portal
             * @tparam object_masks_t the mask container for the objects
             * 
             */
            template <typename object_t, typename object_masks_t>
            class constituents
            {

                friend class volume;

            private:
                vector_type<object_t> _objects;
                object_masks_t _masks;
                transform_store _transforms;

            public:
                /** @return an indexed object with @param object_index - const access */
                const object_t &indexed_object(dindex object_index) const
                {
                    return _objects[object_index];
                }

                /** @return all the objects - const access */
                const vector_type<object_t> &objects() const { return _objects; }

                /** @return the object masks - const access */
                const object_masks_t &masks() const { return _masks; }

                /** @return the object transforms = const access */
                const transform_store &transforms() const { return _transforms; }
            };

            /** Deleted constructor */
            volume() = delete;

            /** Allowed constructors
             * @param name of the volume 
             * 
             * @note will be contructed boundless
             */
            volume(const std::string &name) : _name(name){};
            /** Contructor with name and bounds 
             * @param name of the volume
             * @param bounds of the volume
             */
            volume(const std::string &name, const array_type<scalar, 6> &bounds) : _name(name), _bounds(bounds){};
            volume(const volume &) = default;

            /** @return the bounds - const access */
            const array_type<scalar, 6> &bounds() const { return _bounds; }

            /** @return the name */
            const std::string &name() const { return _name; }

            /** @return the index */
            dindex index() const { return _index; }

            /** @return if the volume is empty or not */
            bool empty() const { return _surfaces.objects().empty(); }

            /** Add a new full set of alignable transforms for surfaces - move semantics
             *
             * @param ctx The context of the call
             * @param trfs The transform container, move semantics
             *
             * @note can throw an exception if input data is inconsistent
             */
            void add_surface_transforms(
                const typename alignable_store::context &ctx,
                typename alignable_store::storage &&trfs) noexcept(false)
            {
                _surfaces._transforms.add_contextual_transforms(ctx, std::move(trfs));
            }

            /** Add the surfaces and their masks - move semantics
             *
             * @param surfaces The (complete) volume surfaces
             * @param surface_mask_container The (complete) surface masks
             */
            void add_surface_components(surface_container &&volume_surfaces, surface_mask_container &&volume_surface_masks)
            {
                _surfaces._objects = std::move(volume_surfaces);
                _surfaces._masks = std::move(volume_surface_masks);
            }

            /** @return all surfaces - const access */
            const auto &surfaces() const { return _surfaces; }

            /** Add the portals, their transforms and their masks, move semantics
             *
             * @param portals The volume (complete set of) portals
             * @param portal_masks The (complete set of) portal masks
             */
            void add_portal_components(portal_container &&portals,
                                       portal_mask_container &&portal_masks)
            {
                _portals._objects = std::move(portals);
                _portals._masks = std::move(portal_masks);
            }

            /** Add a new full set of alignable transforms for portals - move semantics
             *
             * @param ctx The context of the call
             * @param portal_transforms The transform container, move semantics
             *
             * @note can throw an exception if input data is inconsistent
             */
            void add_portal_transforms(
                const typename alignable_store::context &ctx,
                typename alignable_store::storage &&portal_transforms) noexcept(false)
            {
                _portals._transforms.add_contextual_transforms(ctx, std::move(portal_transforms));
            }

            /** @return all portals - const access */
            const auto &portals() const { return _portals; }

        private:
            /** Volume section: name */
            std::string _name = "unknown";
            /** Volume index */
            dindex _index = dindex_invalid;

            /** Bounds section, default for r, z, phi */
            array_type<scalar, 6> _bounds = {0.,
                                             std::numeric_limits<scalar>::max(),
                                             -std::numeric_limits<scalar>::max(),
                                             std::numeric_limits<scalar>::max(),
                                             -M_PI, M_PI};

            /** Surface section */
            constituents<surface, surface_mask_container> _surfaces;

            /** Portal section */
            constituents<portal, portal_mask_container> _portals;
        };

        /** Allowed costructor
         * @param name the detector
         */
        detector(const std::string &name) : _name(name) {}
        detector(const detector & /*ignored*/) = default;
        detector() = delete;
        ~detector() = default;

        /** Add a new volume and retrieve a reference to it
         *
         * @param name of the volume
         * @param bounds of the volume
         *
         * @return non-const reference of the new volume
         */
        volume &new_volume(const std::string &name, const array_type<scalar, 6> &bounds)
        {
            _volumes.push_back(std::move(volume(name, bounds)));
            dindex cvolume_idx = _volumes.size() - 1;
            volume &cvolume = _volumes[cvolume_idx];
            cvolume._index = cvolume_idx;
            return cvolume;
        }

        /** @return the name of the detector */
        const std::string &name() const { return _name; }

        /** @return the contained volumes of the detector - const access */
        const vector_type<volume> &volumes() const { return _volumes; }

        /** @return the volume by @param volume_index - const access */
        const volume &indexed_volume(dindex volume_index) const { return _volumes[volume_index]; }

        /** @return the volume by @param position - const access */
        const volume &indexed_volume(const point3 &p) const
        {
            point2 p2 = {getter::perp(p), p[2]};
            dindex volume_index = _volume_grid.bin(p2);
            return _volumes[volume_index];
        }

        /** @return the volume by @param volume_index - non-const access */
        volume &indexed_volume(dindex volume_index) { return _volumes[volume_index]; }

        /** Add the volume grid - move semantics 
         * 
         * @param v_grid the volume grid to be added
         */
        void add_volume_grid(volume_grid &&v_grid)
        {
            _volume_grid = std::move(v_grid);
        }

        /** @return the volume grid - const access */
        const volume_grid &volume_search_grid() const { return _volume_grid; }

        /** Add local surface finders linked to from the portals - move semantics
         * 
         * This connects portals and surface grids
         */
        void add_surface_finders(vector_type<surface_finder> &&surface_finders)
        {
            _surface_finders = std::move(surface_finders);
        }

        /** @return the surface finders - const access */
        const vector_type<surface_finder>& surface_finders() const { return _surface_finders; }

        /** Output to string */
        const std::string to_string() const
        {
            std::stringstream ss;
            ss << "[>] Detector '" << _name << "' has " << _volumes.size() << " volumes." << std::endl;
            for (const auto &[i, v] : enumerate(_volumes))
            {
                ss << "[>>] Volume at index " << i << " - name: '" << v.name() << "'" << std::endl;
                ss << "     contains    " << v._surfaces.objects().size() << " detector surfaces" << std::endl;
                ss << "                 " << v._portals.objects().size() << " detector portals" << std::endl;
                ss << "     bounds r = (" << v._bounds[0] << ", " << v._bounds[1] << ")" << std::endl;
                ss << "            z = (" << v._bounds[2] << ", " << v._bounds[3] << ")" << std::endl;
            }
            return ss.str();
        };

    private:
        std::string _name = "unknown_detector";
        vector_type<volume> _volumes = {};

        vector_type<surface_finder> _surface_finders;

        volume_grid _volume_grid = volume_grid(std::move(axis::irregular{{}}), std::move(axis::irregular{{}}));
    };

}
