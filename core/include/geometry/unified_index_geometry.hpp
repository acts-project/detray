/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <utility>

#include "core/mask_store.hpp"
#include "geometry/surface_base.hpp"
#include "geometry/volume.hpp"
#include "masks/masks.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {
/**
 * @brief Indexed geometry implementation that unifies surface and portal types
 *
 * This class provides a geometry that defines logic volumes which contain
 * the detector surfaces, joined together by conceptual portal surfaces. It
 * exports all types needed for navigation and strictly only keeps the
 * index data (links) that define the geometry relations. This geometry
 * implemenatation makes no distinction between surfaces and portals. Both
 * carry the same link type: a portal points to the next volume, a surface to
 * the current volume.
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
class unified_index_geometry {

    public:
    // Known primitives
    enum known_objects : bool {
        e_surface = true,
        e_portal = false,
        e_any = false,  // defaults to portal
    };

    /** Encodes the position in a collection container for the respective
        mask type . */
    enum mask_id : unsigned int {
        e_mask_types = 5,
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 3,  // no distinction from surface cylinder
        e_portal_ring2 = 4,
        e_single3 = std::numeric_limits<unsigned int>::max(),
        e_unknown = std::numeric_limits<unsigned int>::max(),
    };

    // Volume type
    using volume_type = volume<array_type>;

    /// volume index: volume the surface belongs to
    using volume_index = dindex;
    /// transform link: transform entry belonging to surface
    using transform_link = dindex;
    /// mask index: type, range
    using mask_index = array_type<dindex, 2>;
    /// edge links: next volume, next (local) object finder
    using edge_links = array_type<dindex, 2>;
    // edge_links are the same for portal and surface and not stored in the
    // masks
    using surface_links = edge_links;
    using portal_links = edge_links;
    /// source link
    using source_link = surface_source_link;
    using bounds_link = source_link;

    /// mask types
    using rectangle = rectangle2<planar_intersector, __plugin::cartesian2,
                                 portal_links, e_rectangle2>;
    using trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, void,
                                 e_trapezoid2>;
    using annulus =
        annulus2<planar_intersector, __plugin::cartesian2, void, e_annulus2>;
    using cylinder = cylinder3<false, cylinder_intersector,
                               __plugin::cylindrical2, void, e_cylinder3>;
    using disc = ring2<planar_intersector, __plugin::cartesian2, void, e_ring2>;

    using mask_container = mask_store<tuple_type, vector_type, rectangle,
                                      trapezoid, annulus, cylinder, disc>;

    /** The Surface definition:
     *  <transform_link, mask_link, volume_link, source_link, edge_link>
     */
    using surface = surface_base<transform_link, mask_index, volume_index,
                                 source_link, edge_links>;
    using surface_container = vector_type<surface>;
    // No difference between surfaces and portals
    using portal = surface;
    using portal_container = surface_container;

    /** Temporary container structures that are used to fill the geometry.
     * The respective objects are sorted by mask type, so that they can be
     * unrolled and filled in lockstep with the masks
     */
    using surface_filling_container =
        array_type<vector_type<surface>, e_mask_types>;
    using portal_filling_container = surface_filling_container;

    /** Default constructor */
    unified_index_geometry() = default;

    /** Copy constructor
     *
     * @param other unified_index_geometry to be copied
     */
    // unified_index_geometry(const unified_index_geometry &other) = default;

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
    template <bool get_surface = true>
    inline size_t n_objects() const {
        return _objects.size();
    }

    /** @return all surfaces/portals in the geometry */
    template <bool get_surface = true>
    inline const auto &objects() const {
        return _objects;
    }

    /** Update the mask links of an object when filling into a large container
     *
     * @param obj the surface or portal
     * @param mask_offset the offset that will be added to the mask links
     */
    inline void update_mask_link(surface &obj, const dindex offset) {
        std::get<1>(obj.mask()) += offset;
    }

    /** Update the transform link of an objects when filling into a large
     * container
     *
     * @param obj the surface or portal
     * @param trsf_offset the offset that will be added to the links
     */
    inline void update_transform_link(surface &obj, const dindex offset) {
        obj.transform() += offset;
    }

    /** Add objects (surfaces/portals) to the geometry
     *
     * @param volume the volume the objects belong to
     * @param surfaces the surfaces that will be filled into the volume
     */
    template <bool add_surfaces = true>
    inline void add_objects(volume_type &volume,
                            const surface_container &surfaces) {
        const auto offset = _objects.size();
        _objects.reserve(_objects.size() + surfaces.size());
        _objects.insert(_objects.end(), surfaces.begin(), surfaces.end());

        volume.template set_range<add_surfaces>({offset, _objects.size()});
    }

    private:
    /** Contains the geometrical relations*/
    vector_type<volume_type> _volumes = {};

    /** All surfaces and portals in the geometry in contigous memory */
    surface_container _objects = {};
};

}  // namespace detray