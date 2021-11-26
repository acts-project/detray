/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <map>
#include <sstream>
#include <string>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/core/mask_store.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/object_registry.hpp"
#include "detray/geometry/surface_base.hpp"
#include "detray/geometry/volume.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tools/local_object_finder.hpp"

namespace detray {

// Algebra, point2 is not strongly typed
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

/** The detector definition.
 *
 * This class is a heavy templated detector definition class, that sets the
 * interface between geometry, navigator and grid.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 * @tparam tuple_type the type of the internal tuple, must have STL semantics
 * @tparam vector_type the type of the internal array, must have STL semantics
 * @tparam alignable_store the type of the transform store
 * @tparam geometry_type the geometry implementation to be used
 * @tparam surfaces_populator_type the type of populator used to fill the
 * surfaces grids
 * @tparam surfaces_serializer_type the type of the memory serializer for the
 * surfaces grids
 *
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          template <typename...> class jagged_vector_type = djagged_vector,
          typename surfaces_serializer_type = serializer2,
          typename name_map = std::map<dindex, std::string>,
          typename source_link = dindex>
class detector {

    public:
    /// Forward the alignable container and context
    using transform_store = static_transform_store<vector_type>;
    using context = typename transform_store::context;

    // TODO: Remove this from detector
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

    /// volume index: volume the surface belongs to
    using volume_link = dindex;
    /// transform link: transform entry belonging to surface
    using transform_link = dindex;
    /// mask index: type, range
    using mask_link = array_type<dindex, 2>;
    /// edge links: next volume, next (local) object finder
    using edge_type = array_type<dindex, 2>;

    /// mask types
    using rectangle =
        rectangle2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_type, e_rectangle2>;
    using trapezoid =
        trapezoid2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_type, e_trapezoid2>;
    using annulus =
        annulus2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                 edge_type, e_annulus2>;
    using cylinder = cylinder3<false, cylinder_intersector,
                               __plugin::cylindrical2<detray::scalar>,
                               edge_type, e_cylinder3>;
    using disc = ring2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                       edge_type, e_ring2>;

    using mask_container = mask_store<tuple_type, vector_type, rectangle,
                                      trapezoid, annulus, cylinder, disc>;

    /** The Surface definition:
     *  <transform_link, mask_link, volume_link, source_link, edge_link>
     */
    using surface_type = surface_base<transform_link, mask_link, volume_link,
                                      source_link, edge_type>;

    using object_id = object_registry<surface_type>;
    using surface_container = vector_type<surface_type>;

    // Volume type
    using volume_type = volume<object_id, dindex_range, array_type>;

    /// Accelerator structure

    /// Volume grid definition
    using volume_grid =
        grid2<replace_populator, axis::irregular, axis::irregular, serializer2,
              vector_type, jagged_vector_type, array_type, tuple_type, dindex>;

    using surfaces_regular_circular_grid =
        grid2<attach_populator, axis::regular, axis::circular,
              surfaces_serializer_type, vector_type, jagged_vector_type,
              array_type, tuple_type, dindex, false>;

    using surfaces_regular_axis =
        typename surfaces_regular_circular_grid::axis_p0_t;
    using surfaces_circular_axis =
        typename surfaces_regular_circular_grid::axis_p1_t;

    // Neighborhood finder, using accelerator data structure
    using surfaces_finder = surfaces_regular_circular_grid;

    /** Temporary container structures that are used to fill the detector.
     * The respective objects are sorted by mask type, so that they can be
     * unrolled and filled in lockstep with the masks
     */
    using surface_filling_container =
        array_type<vector_type<surface_type>, e_mask_types>;
    using transform_container =
        array_type<transform_store, mask_id::e_mask_types>;

    detector() = delete;

    /** Allowed costructor
     * @param resource memory resource for the allocation of members
     */
    DETRAY_HOST
    detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _volume_grid(std::move(typename volume_grid::axis_p0_t{resource}),
                       std::move(typename volume_grid::axis_p1_t{resource}),
                       resource),
          //_surfaces_finders(&resource),
          _resource(&resource) {}

    /** Constructor with detector_data **/
    template <typename detector_data_type,
              std::enable_if_t<!std::is_base_of_v<vecmem::memory_resource,
                                                  detector_data_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE detector(detector_data_type &det_data)
        : _volumes(det_data._volumes_data),
          _surfaces(det_data._surfaces_data),
          _transforms(det_data._transforms_data),
          _masks(det_data._masks_data),
          _volume_grid(det_data._volume_grid_view) {}

    /** Add a new volume and retrieve a reference to it
     *
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    DETRAY_HOST
    volume_type &new_volume(const array_type<scalar, 6> &bounds,
                            dindex surfaces_finder_entry = dindex_invalid) {
        volume_type &cvolume = _volumes.emplace_back(bounds);
        cvolume.set_index(_volumes.size() - 1);
        cvolume.set_surfaces_finder(surfaces_finder_entry);

        return cvolume;
    }

    /** @return the contained volumes of the detector - const access */
    DETRAY_HOST_DEVICE
    inline auto &volumes() const { return _volumes; }

    /** @return the contained volumes of the detector - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &volumes() { return _volumes; }

    /** @return the volume by @param volume_index - const access */
    DETRAY_HOST_DEVICE
    inline auto &volume_by_index(dindex volume_index) const {
        return _volumes[volume_index];
    }

    /** @return the volume by @param volume_index - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &volume_by_index(dindex volume_index) {
        return _volumes[volume_index];
    }

    /** @return the volume by @param position - const access */
    DETRAY_HOST_DEVICE
    inline auto &volume_by_pos(const point3 &p) const {
        point2 p2 = {getter::perp(p), p[2]};
        dindex volume_index = _volume_grid.bin(p2);
        return _volumes[volume_index];
    }

    /** @return all objects of a given type - const access */
    DETRAY_HOST_DEVICE
    inline auto &surfaces() const { return _surfaces; }

    /** @return all objects of a given type - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &surfaces() { return _surfaces; }

    /** @return all surface/portal masks in the geometry - const access */
    DETRAY_HOST_DEVICE
    inline auto &masks() const { return _masks; }

    /** @return all surface/portal masks in the geometry - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &masks() { return _masks; }

    /** Add pre-built mask store
     *
     * @param masks the conatiner for surface masks
     */
    DETRAY_HOST
    inline void add_mask_store(mask_container &&masks) {
        _masks = std::move(masks);
    }

    /** Get all transform in an index range from the detector
     *
     * @param range The range of surfaces/portals in the transform store
     * @param ctx The context of the call
     *
     * @return ranged iterator to the object transforms
     */
    DETRAY_HOST_DEVICE
    inline const auto transforms(const dindex_range &range,
                                 const context &ctx = {}) const {
        return _transforms.range(range, ctx);
    }

    /** Get all transform in an index range from the detector - const
     *
     * @param ctx The context of the call
     *
     * @return detector transform store
     */
    DETRAY_HOST_DEVICE
    inline const auto &transforms(const context &ctx = {}) const {
        return _transforms;
    }

    DETRAY_HOST_DEVICE
    inline auto &transforms(const context &ctx = {}) { return _transforms; }

    /** Add pre-built transform store
     *
     * @param transf the constianer for surface transforms
     */
    DETRAY_HOST
    inline void add_transform_store(transform_store &&transf) {
        _transforms = std::move(transf);
    }

    /** Get all available data from the detector without std::tie
     *
     * @param ctx The context of the call
     *
     * @return a struct that contains references to all relevant containers.
     */
    DETRAY_HOST_DEVICE
    const auto data(const context &ctx = {}) const {
        struct data_core {
            const dvector<volume_type> &volumes;
            const transform_store &transforms;
            const mask_container &masks;
            const surface_container &surfaces;
        };
        return data_core{_volumes, _transforms, _masks, _surfaces};
    }

    /** Add a new full set of detector components (e.g. transforms or volumes)
     *  according to given context.
     *
     * @tparam detector_components types of detector components
     * @tparam object_type check whether we deal with e.g. surfaces or portals
     *
     * @param ctx The context of the call
     * @param components The components to be added
     *
     * @note can throw an exception if input data is inconsistent
     */
    template <typename... detector_components>
    DETRAY_HOST inline void add_objects(
        const context ctx,
        detector_components &&... components) noexcept(false) {
        // Fill according to type, starting at type '0' (see 'mask_id')
        fill_containers(ctx, std::forward<detector_components>(components)...);
    }

    /** Unrolls the data containers according to the mask type and fill the
     *  global containers. It registers the indexing in the geometry.
     *
     * @tparam current_type the current mask context to be processed
     * @tparam surface_container surfaces/portals for which the links are
     * updated
     * @tparam mask_container surface/portal masks, sorted by type
     * @tparam object_type check whether we deal with surfaces or portals
     *
     * @param volume The volume we add the transforms to
     * @param objects The geometry objects in the volume
     * @param trfs The transforms
     * @param ctx The context of the call
     *
     * @note can throw an exception if input data is inconsistent
     */
    template <unsigned int current_type = 0, typename surface_container>
    DETRAY_HOST inline void fill_containers(
        const context ctx, volume_type &volume, surface_container &surfaces,
        mask_container &masks, transform_container &trfs) noexcept(false) {
        // Get the surfaces/portals for a mask type
        auto &typed_surfaces = surfaces[current_type];
        // Get the corresponding transforms
        const auto &object_transforms = trfs[current_type];
        // and the corresponding masks
        auto &object_masks = masks.template group<current_type>();

        if (not object_transforms.empty(ctx) and not typed_surfaces.empty()) {
            // Current offsets into detectors containers
            const auto trsf_offset = _transforms.size(ctx);
            const auto mask_offset = _masks.template size<current_type>();

            // Fill the correct mask type
            _masks.add_masks(object_masks);
            _transforms.append(ctx, std::move(std::get<current_type>(trfs)));

            // Update the surfaces mask link
            for (auto &obj : typed_surfaces) {
                detail::get<1>(obj.mask()) += mask_offset;
                obj.transform() += trsf_offset;
            }

            // Now put the updated objects into the geometry
            const auto offset = _surfaces.size();
            _surfaces.reserve(_surfaces.size() + typed_surfaces.size());
            _surfaces.insert(_surfaces.end(), typed_surfaces.begin(),
                             typed_surfaces.end());

            volume.set_range({offset, _surfaces.size()});
        }

        // Next mask type
        if constexpr (current_type <
                      std::tuple_size_v<typename mask_container::mask_tuple> -
                          1) {
            return fill_containers<current_type + 1, surface_container>(
                ctx, volume, surfaces, masks, trfs);
        }
        // If no mask type fits, don't fill the data.
    }

    /** Add the volume grid - move semantics
     *
     * @param v_grid the volume grid to be added
     */
    DETRAY_HOST
    inline void add_volume_grid(volume_grid &&v_grid) {
        _volume_grid = std::move(v_grid);
    }

    /** @return the volume grid - const access */
    DETRAY_HOST_DEVICE
    inline const volume_grid &volume_search_grid() const {
        return _volume_grid;
    }

    DETRAY_HOST_DEVICE
    inline volume_grid &volume_search_grid() { return _volume_grid; }

    /** Add local surface finders linked to from the portals - move semantics
     *
     * This connects portals and surface grids
     */
    DETRAY_HOST
    inline void add_surfaces_finders(
        vector_type<surfaces_finder> &&surfaces_finders) {
        //_surfaces_finders = std::move(surfaces_finders);
    }

    /** @return the surface finders - const access */
    DETRAY_HOST_DEVICE
    inline const vector_type<surfaces_finder> &surfaces_finders() const {
        // return _surfaces_finders;
    }

    /** Output to string */
    DETRAY_HOST
    const std::string to_string(const name_map &names) const {
        std::stringstream ss;

        ss << "[>] Detector '" << names.at(0) << "' has " << _volumes.size()
           << " volumes." << std::endl;
        // ss << "    contains  " << _surfaces_finders.size()
        //   << " local surface finders." << std::endl;

        for (const auto &[i, v] : enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << ": " << std::endl;
            ss << " - name: '" << v.name(names) << "'" << std::endl;

            ss << "     contains    "
               << v.template n_objects<object_id::e_surface>() << " surfaces "
               << std::endl;

            ss << "                 "
               << v.template n_objects<object_id::e_portal>() << " portals "
               << std::endl;

            if (v.surfaces_finder_entry() != dindex_invalid) {
                ss << "  sf finders idx " << v.surfaces_finder_entry()
                   << std::endl;
            }

            const auto &bounds = v.bounds();
            ss << "     bounds r = (" << bounds[0] << ", " << bounds[1] << ")"
               << std::endl;
            ss << "            z = (" << bounds[2] << ", " << bounds[3] << ")"
               << std::endl;
        }

        return ss.str();
    }

    /** @return the pointer of memoery resource */
    DETRAY_HOST
    auto resource() { return _resource; }

    private:
    /** Contains the geometrical relations */
    vector_type<volume_type> _volumes;

    /** All surfaces and portals in the geometry in contiguous memory */
    surface_container _surfaces;

    /** Keeps all of the transform data in contiguous memory*/
    transform_store _transforms;

    /** Surface and portal masks of the detector in contiguous memory */
    mask_container _masks;

    /* TODO: surfaces_finder needs to be refactored */
    // vecmem::vector<surfaces_finder> _surfaces_finders;

    volume_grid _volume_grid;

    vecmem::memory_resource *_resource = nullptr;
};

/** A static inplementation of detector data for device
 *
 */
template <typename detector_type>
struct detector_data {

    // type definitions
    using volume_t = typename detector_type::volume_type;
    using surface_t = typename detector_type::surface_type;
    using mask_container_t = typename detector_type::mask_container;
    using transform_store_t = typename detector_type::transform_store;
    using volume_grid_t = typename detector_type::volume_grid;

    detector_data(detector_type &det)
        : _volumes_data(vecmem::get_data(det.volumes())),
          _surfaces_data(vecmem::get_data(det.surfaces())),
          _masks_data(get_data(det.masks())),
          _transforms_data(get_data(det.transforms())),
          _volume_grid_data(
              get_data(det.volume_search_grid(), *det.resource())),
          _volume_grid_view(_volume_grid_data) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    mask_store_data<mask_container_t> _masks_data;
    static_transform_store_data<transform_store_t> _transforms_data;
    grid2_data<volume_grid_t> _volume_grid_data;
    grid2_view<volume_grid_t> _volume_grid_view;
};

/** stand alone function for detector_data get function
 **/
template <template <typename, unsigned int> class array_type,
          template <typename...> class tuple_type,
          template <typename...> class vector_type,
          template <typename...> class jagged_vector_type,
          typename surfaces_serializer_type, typename name_map,
          typename source_link>
inline detector_data<
    detector<array_type, tuple_type, vector_type, jagged_vector_type,
             surfaces_serializer_type, name_map, source_link> >
get_data(detector<array_type, tuple_type, vector_type, jagged_vector_type,
                  surfaces_serializer_type, name_map, source_link> &det) {
    return det;
}

}  // namespace detray
