/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <type_traits>
#include <utility>

#include "core/mask_store.hpp"
#include "geometry/surface_base.hpp"
#include "geometry/volume.hpp"
#include "masks/masks.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {
/**
 * @brief Index geometry implementation
 *
 * This class provides a geometry that defines logic volumes which contain
 * the detector surfaces, joined together by dedicated portal surfaces. It
 * exports all types needed for navigation and strictly only keeps the
 * index data (links) that define the geometry relations.
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
template <template <typename, unsigned int> class array_type = darray,
          template <typename> class vector_type = dvector,
          template <typename...> class tuple_type = dtuple,
          typename surface_source_link = dindex,
          typename bounds_source_link = dindex>
class index_geometry {

    public:
    // Known primitives
    enum object_id : unsigned int {
        e_object_types = 2,
        e_surface = 0,
        e_portal = 1,
        e_any = 1,  // defaults to portal
    };

    /** Encodes the position in a collection container for the respective
        mask type . */
    enum mask_id : unsigned int {
        e_mask_types = 6,
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 4,
        e_portal_ring2 = 5,
        e_single3 = std::numeric_limits<unsigned int>::max(),
        e_unknown = std::numeric_limits<unsigned int>::max(),
    };

    // Volume type
    using volume_type = volume<object_id, dindex_range, array_type>;

    /// Portals components:
    /// - links:  next volume, next (local) object finder
    using portal_links = array_type<dindex, 2>;
    /// - masks, with mask identifiers 0, 1
    using portal_cylinder =
        cylinder3<false, cylinder_intersector, __plugin::cylindrical2,
                  portal_links, e_portal_cylinder3>;
    using portal_disc = ring2<planar_intersector, __plugin::cartesian2,
                              portal_links, e_portal_ring2>;
    // - mask index: type, { first/last }
    using portal_mask_index = tuple_type<dindex, array_type<dindex, 2>>;

    /** The Portal definition:
     *  <transform_link, mask_index, volume_link, source_link >
     *
     * transform_link: index into the transform container
     * mask_index: typed index into the mask container
     * volume_link: index of the volume this portal belongs to
     * source_link: some link to an eventual exernal representation
     *
     */
    using bounds_link = bounds_source_link;
    using portal = surface_base<dindex, portal_mask_index, dindex, bounds_link,
                                portal_links>;
    using portal_container = vector_type<portal>;

    /// Surface components:
    /// - surface links
    using surface_links = array_type<dindex, 1>;
    /// - masks, with mask identifiers 0,1,2
    using surface_rectangle =
        rectangle2<planar_intersector, __plugin::cartesian2, surface_links,
                   e_rectangle2>;
    using surface_trapezoid =
        trapezoid2<planar_intersector, __plugin::cartesian2, surface_links,
                   e_trapezoid2>;
    using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2,
                                     surface_links, e_annulus2>;
    using surface_cylinder =
        cylinder3<false, cylinder_intersector, __plugin::cylindrical2,
                  surface_links, e_cylinder3>;
    /// - mask index: type, entry
    using surface_mask_index = array_type<dindex, 2>;

    using mask_container =
        mask_store<tuple_type, vector_type, surface_rectangle,
                   surface_trapezoid, surface_annulus, surface_cylinder,
                   portal_cylinder, portal_disc>;

    using source_link = surface_source_link;
    /** The Surface definition:
     *  <transform_link, mask_link, volume_link, source_link >
     */
    using surface = surface_base<dindex, surface_mask_index, dindex,
                                 source_link, surface_links>;
    using surface_container = vector_type<surface>;

    /** Temporary container structures that are used to fill the geometry.
     * The respective objects are sorted by mask type, so that they can be
     * unrolled and filled in lockstep with the masks
     */
    using surface_filling_container =
        array_type<vector_type<surface>, e_mask_types>;
    using portal_filling_container =
        array_type<vector_type<portal>, e_mask_types>;

    /** Default constructor */
    index_geometry() = default;

    /** @return total number of volumes */
    const size_t n_volumes() const { return _volumes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto &volumes() const { return _volumes; }

    /** @return the volume by @param volume_index - const access. */
    inline const volume_type &volume_by_index(dindex volume_index) const {
        return _volumes[volume_index];
    }

    /** @return the volume by @param volume_index - non-const access. */
    inline volume_type &volume_by_index(dindex volume_index) {
        return _volumes[volume_index];
    }

    /** Add a new volume and retrieve a reference to it
     *
     * @param name of the volume
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    inline volume_type &new_volume(
        const array_type<scalar, 6> &bounds,
        dindex surfaces_finder_entry = dindex_invalid) {
        volume_type &cvolume = _volumes.emplace_back(bounds);
        cvolume.set_index(_volumes.size() - 1);
        cvolume.set_surfaces_finder(surfaces_finder_entry);

        return cvolume;
    }

    /** @return all surfaces/portals in the geometry */
    template <object_id object_type = e_surface>
    inline size_t n_objects() const {
        if constexpr (object_type == e_surface) {
            return _surfaces.size();
        } else {
            return _portals.size();
        }
    }

    /** @return all surfaces/portals in the geometry */
    template <object_id object_type = e_surface>
    inline constexpr const auto &objects() const {
        if constexpr (object_type == e_surface) {
            return _surfaces;
        } else {
            return _portals;
        }
    }

    /** Update the mask link of a surface when filling into a large container
     *
     * @param sf the surface
     * @param mask_offset the offset that will be added to the mask links
     */
    inline void update_mask_link(surface &sf, const dindex mask_offset) {
        std::get<1>(sf.mask()) += mask_offset;
    }

    /** Update the mask links of a portal when filling into a large container
     *
     * @param pt the portal
     * @param mask_offset the offset that will be added to the mask links
     */
    inline void update_mask_link(portal &pt, const dindex mask_offset) {
        auto &portal_mask_index = std::get<1>(pt.mask());
        portal_mask_index[0] += mask_offset;
        portal_mask_index[1] += mask_offset;
    }

    /** Update the transform link of a surface/portal when filling into a large
     * container
     *
     * @param sf the surface
     * @param trsf_offset the offset that will be added to the links
     */
    template <typename object_t>
    inline void update_transform_link(object_t &obj, const dindex trsf_offset) {
        obj.transform() += trsf_offset;
    }

    /** Add objects (surfaces/portals) to the geometry
     *
     * @param volume the volume the objects belong to
     * @param surfaces the surfaces that will be filled into the volume
     */
    template <typename object_t>
    inline void add_objects(volume_type &volume,
                            const vector_type<object_t> &objects) {
        if constexpr (std::is_same_v<object_t, surface>) {
            const auto offset = _surfaces.size();
            _surfaces.reserve(_surfaces.size() + objects.size());
            _surfaces.insert(_surfaces.end(), objects.begin(), objects.end());

            volume.template set_range<e_surface>({offset, _surfaces.size()});
        } else {
            const auto offset = _portals.size();
            _portals.reserve(_portals.size() + objects.size());
            _portals.insert(_portals.end(), objects.begin(), objects.end());

            volume.template set_range<e_portal>({offset, _portals.size()});
        }
    }

    private:
    /** Contains the geometrical relations */
    vector_type<volume_type> _volumes = {};

    /** All surfaces and portals in the geometry in contigous memory */
    surface_container _surfaces = {};
    portal_container _portals = {};
};

}  // namespace detray