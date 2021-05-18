/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "core/surface.hpp"
#include "core/transform_store.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "masks/annulus2.hpp"
#include "masks/rectangle2.hpp"
#include "masks/trapezoid2.hpp"
#include "masks/cylinder3.hpp"
#include "masks/ring2.hpp"
#include "utils/indexing.hpp"
#include "utils/containers.hpp"
#include "utils/enumerate.hpp"
#include "tools/planar_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "tools/concentric_cylinder_intersector.hpp"

#include <string>
#include <sstream>

namespace detray
{

    // Algebra, point2 is not strongly typed
    using point3 = __plugin::transform3::point3;
    using vector3 = __plugin::transform3::vector3;
    using point2 = __plugin::cartesian2::point2;

    /** Indexed detector definition.
     *
     * @tparam alignable_store the type of the transform store
     * @tparam surface_source_link the type of the link to an external surface source
     * @tparam bounds_source_link the type of the link to an external bounds source
     * 
     */
    template <typename alignable_store = static_transform_store,
              typename surface_source_link = dindex,
              typename bounds_source_link = dindex>
    class proto_detector
    {

    public:
        /// Volume grid definition
        using volume_grid = grid2<replace_populator<>, axis::irregular, axis::irregular, serializer2>;

        /// Portals components:
        /// - links:  next volume, next object finder
        using portal_links = darray<dindex, 2>;
        /// - masks, with mask identifiers 0, 1
        using portal_cylinder = cylinder3<false, concentric_cylinder_intersector, __plugin::cylindrical2, portal_links, 0>;
        using portal_disc = ring2<planar_intersector, __plugin::cartesian2, portal_links, 1>;
        // - mask index: type, first/last
        using portal_mask_index = darray<dindex, 3>;
        using portal_mask_container = dtuple<dvector<portal_cylinder>, dvector<portal_disc>>;
        // - portal transforms, different from surfaces w/o alignment
        using portal_transforms = dvector<transform3>;

        /// The Portal definition:
        ///  <transform_link, mask_index, volume_link, source_link >
        using portal = surface<dindex, portal_mask_index, dindex, surface_source_link>;
        using portals = dvector<portal>;

        /// Surface components:
        /// - masks, with mask identifiers 0,1,2
        using surface_rectangle = rectangle2<planar_intersector, __plugin::cartesian2, bounds_source_link, 0>;
        using surface_trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, bounds_source_link, 1>;
        using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2, bounds_source_link, 2>;
        using surface_cylinder = cylinder3<true, cylinder_intersector, __plugin::cylindrical2, bounds_source_link, 3>;
        /// - mask index: type, entry
        using surface_mask_index = darray<dindex, 2>;
        using surface_mask_container = dtuple<dvector<surface_rectangle>,
                                     dvector<surface_trapezoid>,
                                     dvector<surface_annulus>,
                                     dvector<surface_cylinder>>;

        /** The Surface definition:
         *  <transform_link, mask_link, volume_link, source_link >
         */
        using surface = surface<dindex, surface_mask_index, dindex, surface_source_link>;
        using surface_container = dvector<surface>;

        /** Nested volume struct that holds the local information of the
         * volume and its portals.
         */
        class volume
        {

            friend class proto_detector<alignable_store, surface_source_link, bounds_source_link>;

        public:
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
            volume(const std::string &name, const darray<scalar, 6> &bounds) : _name(name), _bounds(bounds){};
            volume(const volume &) = default;

            /** @return the bounds - const access */
            const darray<scalar, 6> &bounds() const { return _bounds; }

            /** @return the name */
            const std::string &name() const { return _name; }

            /** @return the index */
            dindex index() const { return _index; }

            /** @return if the volume is empty or not */
            bool empty() const { return _surfaces.empty(); }

            /** Add a new full set of alignable transforms - move semantics
             *
             * @param ctx The context of the call
             * @param trfs The transform container, move semantics
             *
             * @note can throw an exception if input data is inconsistent
             */
            void add_contextual_transforms(
                const typename alignable_store::context &ctx,
                typename alignable_store::storage &&trfs) noexcept(false)
            {
                _surface_transforms.add_contextual_transforms(ctx, std::move(trfs));
            }

            /** Add a new full set of alignable transforms - copy semantics
             *
             * @param ctx The context of the call
             * @param trfs The transform container, move semantics
             *
             * @note can throw an exception if input data is inconsistent
             */
            void add_contextual_transforms(
                const typename alignable_store::context &ctx,
                const typename alignable_store::storage &trfs) noexcept(false)
            {
                _surface_transforms.add_contextual_transforms(ctx, trfs);
            }

            /** Add the surfaces and their masks - move semantics
             *
             * @param surfaces The (complete) volume surfaces
             * @param surface_mask_container The (complete) surface masks
             */
            void add_surface_components(surface_container &&volume_surfaces, surface_mask_container &&volume_surface_mask_container)
            {
                _surfaces = std::move(volume_surfaces);
                _surface_mask_container = std::move(volume_surface_mask_container);
            }

            /** Add the surfaces and their masks - copy semantics
             *
             * @param surfaces The (complete) volume surfaces
             * @param surface_mask_container The (complete) surface masks
             */
            void add_surface_components(const surface_container &volume_surfaces, const surface_mask_container &volume_surface_mask_container)
            {
                _surfaces = volume_surfaces;
                _surface_mask_container = volume_surface_mask_container;
            }

            /** Add the portals, their transforms and their masks, move semantics
             *
             * @param volume_portals The volume portals
             * @param volume_portal_transforms The (complete) portal transforms
             * @param volume_portal_mask_container The (complete) portal masks
             */
            void add_portal_components(portals &&volume_portals,
                                       portal_transforms &&volume_portal_transforms,
                                       portal_mask_container &&volume_portal_mask_container)
            {
                _portals = std::move(volume_portals);
                _portal_transforms = std::move(volume_portal_transforms);
                _portal_mask_container = std::move(volume_portal_mask_container);
            }

            /** Const Access to the detector surfaces */
            const surface_container& surfaces() const { return _surfaces; }

            /** Const Access to the surface transform store */
            const alignable_store& surface_transforms() const { return _surface_transforms; }

            /** Const Access to the surface masks */
            const surface_mask_container& masks() const { return _surface_mask_container; }

        private:
            /// Volume section: name
            std::string _name = "unknown";
            /// Volume index
            dindex _index = dindex_invalid;

            /// Bounds section, default for r, z, phi
            darray<scalar, 6> _bounds = {0.,
                                         std::numeric_limits<scalar>::max(),
                                         -std::numeric_limits<scalar>::max(),
                                         std::numeric_limits<scalar>::max(),
                                         -M_PI, M_PI};

            /// Surface section
            surface_mask_container _surface_mask_container;
            surface_container _surfaces;
            alignable_store _surface_transforms;

            /// Portal section
            portal_mask_container _portal_mask_container;
            portals _portals;
            portal_transforms _portal_transforms;
        };

        /** Allowed costructor
         * @param name the detector
         */
        proto_detector(const std::string &name) : _name(name) {}
        proto_detector(const proto_detector & /*ignored*/) = default;
        proto_detector() = delete;
        ~proto_detector() = default;

        /** Add a new volume and retrieve a reference to it
         *
         * @param name of the volume
         * @param bounds of the volume
         *
         * @return non-const reference of the new volume
         */
        volume &new_volume(const std::string &name, const darray<scalar, 6> &bounds)
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
        const dvector<volume> &volumes() const { return _volumes; }

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

        /** Output to string */
        const std::string to_string() const
        {
            std::stringstream ss;
            ss << "[>] Detector '" << _name << "' has " << _volumes.size() << " volumes." << std::endl;
            for (const auto &[i, v] : enumerate(_volumes))
            {
                ss << "[>>] Volume at index " << i << " - name: '" << v.name() << "'" << std::endl;
                ss << "     contains " << v._surfaces.size() << " detector surfaces" << std::endl;
                ss << "              " << v._portals.size() << " detector portals" << std::endl;
            }
            return ss.str();
        };

    private:
        std::string _name = "unknown_detector";
        dvector<volume> _volumes = {};

        volume_grid _volume_grid = volume_grid(std::move(axis::irregular{{}}), std::move(axis::irregular{{}}));
    };

}
