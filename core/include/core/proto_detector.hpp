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
#include "masks/annulus2.hpp"
#include "masks/rectangle2.hpp"
#include "masks/trapezoid2.hpp"
#include "masks/cylinder3.hpp"
#include "masks/ring2.hpp"
#include "utils/indexing.hpp"
#include "utils/containers.hpp"
#include "utils/enumerate.hpp"
#include "tools/planar_intersector.hpp"
#include "tools/concentric_cylinder_intersector.hpp"

#include <string>

namespace detray
{

    /// Algebra, point2 is not strongly typed
    using point3 = __plugin::transform3::point3;
    using vector3 = __plugin::transform3::vector3;
    using point2 = __plugin::cartesian2::point2;

    /// Indexed detector definition.
    ///
    /// The surface_source_link is
    template <typename alignable_store = static_transform_store,
              typename surface_source_link = dindex,
              typename bounds_source_link = dindex>
    class proto_detector
    {

    public:
        /// Portals components:
        /// - links:  opposite volume, along volume, opposite object finder, along volume finder
        using portal_links = darray<dindex, 4>;
        /// - masks, with mask identifiers 0, 1
        using portal_cylinder = cylinder3<false, concentric_cylinder_intersector, __plugin::cylindrical2, portal_links, 0>;
        using portal_disc = ring2<planar_intersector, __plugin::cartesian2, portal_links, 1>;
        // - mask index: type, first/last
        using portal_mask_index = darray<dindex, 3>;
        using portal_masks = dtuple<dvector<portal_cylinder>, dvector<portal_disc> >;
        // - portal transforms, different from surfaces w/o alignment
        using portal_transforms = dvector<transform3>;

        /// The Portal definition:
        ///  <transform_link, mask_link, volume_link, source_link >
        using portal_surface = surface<dindex, portal_mask_index, dindex, surface_source_link>;
        using portals = dvector<portal_surface>;

        /// Surface components:
        /// - masks, with mask identifiers 0,1,2
        using surface_rectangle = rectangle2<planar_intersector, __plugin::cartesian2, bounds_source_link, 0>;
        using surface_trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, bounds_source_link, 1>;
        using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2, bounds_source_link, 2>;
        /// - mask index: type, entry
        using surface_mask_index = darray<dindex, 2>;
        using surface_masks = dtuple<dvector<surface_rectangle>, dvector<surface_trapezoid>, dvector<surface_annulus> >;

        /// The Surface definition:
        ///  <transform_link, mask_link, volume_link, source_link >
        using detector_surface = surface<dindex, surface_mask_index, dindex, surface_source_link>;
        using surfaces = dvector<detector_surface>;

        /// Nested volume struct that holds the local information of the
        /// volume and its portals.
        class volume
        {

            friend class proto_detector<alignable_store, surface_source_link, bounds_source_link>;

        public:
            /// Deleted constructor
            volume() = delete;

            /// Allowed constructors
            /// @param name the volume
            volume(const std::string &name) : _name(name) {}
            volume(const volume &) = default;

            /// @return the name
            const std::string &name() const { return _name; }

            /// Add a new full set of alignable transforms
            ///
            /// @param ctx The context of the call
            /// @param trfs The transform container, move semantics
            ///
            /// @note can throw an exception if input data is inconsistent
            void add_contextual_transforms(
                const typename alignable_store::context &ctx,
                typename alignable_store::storage &&trfs) noexcept(false)
            {
                _surface_transforms.add_contextual_transforms(ctx, std::move(trfs));
            }

            /// Add the surfaces and their masks, move semantics
            ///
            /// @param surfaces The (complete) volume surfaces
            /// @param surface_masks The (complete) surface masks
            void add_surface_components(surfaces &&volume_surfaces, surface_masks &&volume_surface_masks)
            {
                _surfaces = std::move(volume_surfaces);
                _surface_masks = std::move(volume_surface_masks);
            }

            /// Add the portals, their transforms and their masks, move semantics
            ///
            /// @param volume_portals The volume portals
            /// @param volume_portal_transforms The (complete) portal transforms
            /// @param volume_portal_masks The (complete) portal masks
            void add_portal_components(portals &&volume_portals,
                                       portal_transforms &&volume_portal_transforms,
                                       portal_masks &&volume_portal_masks)
            {
                _portals = std::move(volume_portals);
                _portal_transforms = std::move(volume_portal_transforms);
                _portal_masks = std::move(volume_portal_masks);
            }

        private:
            // Volume section
            std::string _name = "unknown";
            dindex index = dindex_invalid;

            // Surface section
            surface_masks _surface_masks;
            surfaces _surfaces;
            alignable_store _surface_transforms;

            // Portal section
            portal_masks _portal_masks;
            portals _portals;
            portal_transforms _portal_transforms;
        };

        /// Allowed costructor
        /// @param name the detector
        proto_detector(const std::string &name) : _name(name) {}
        proto_detector(const proto_detector & /*ignored*/) = default;

        /// Delete constructor
        proto_detector() = delete;

        /// Default destructor, non-virtual
        ~proto_detector() = default;

        /// Add a new volume and retrieve a reference to it
        ///
        /// @param name the volume
        volume &new_volume(const std::string &name)
        {
            _volumes.push_back(std::move(volume(name)));
            dindex cvolume_idx = _volumes.size() - 1;
            volume &cvolume = _volumes[cvolume_idx];
            cvolume.index = cvolume_idx;
            return cvolume;
        }

        /// @return the name of the detector
        const std::string &name() const { return _name; }

    private:
        std::string _name = "unknown_detector";
        dvector<volume> _volumes = {};
    };

}
