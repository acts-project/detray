/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <map>
#include <sstream>
#include <string>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/core/intersection.hpp"
#include "detray/core/mask_store.hpp"
#include "detray/core/surfaces_finder.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/volume.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"

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
 * @tparam metadata helper that defines collection and link types centrally
 * @tparam array_type the type of the internal array, must have STL semantics
 * @tparam tuple_type the type of the internal tuple, must have STL semantics
 * @tparam vector_type the type of the internal array, must have STL semantics
 * @tparam source_link the surface source link
 */
template <typename metadata,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          typename source_link = dindex>
class detector {

    public:
    template <typename T>
    using vector_type = vector_t<T>;

    using name_map = std::map<dindex, std::string>;

    /// Forward the alignable container and context
    using transform_container =
        typename metadata::template transform_store<vector_t>;
    using transform_link = typename transform_container::link_type;
    using context = typename transform_container::context;

    /// Forward mask types
    using masks = typename metadata::mask_definitions;
    using mask_container =
        typename masks::template container_type<tuple_t, vector_t>;

    /// volume index: volume the surface belongs to
    using volume_link = dindex;
    using surface_type =
        surface<masks, transform_link, volume_link, source_link>;

    using objects =
        typename metadata::template object_definitions<surface_type>;
    using surface_container = vector_t<surface_type>;
    // Volume type
    using volume_type = volume<objects, dindex_range, array_t>;

    /** Temporary container structures that are used to fill the detector.
     * The respective objects are sorted by mask type, so that they can be
     * unrolled and filled in lockstep with the masks
     */
    // TODO: Move to volume builder in the future
    using surface_filling_container =
        array_t<vector_t<surface_type>, masks::n_types>;
    using transform_filling_container =
        array_t<transform_container, masks::n_types>;

    /// Accelerator structures

    /// Volume finder definition
    using volume_finder =
        typename metadata::template volume_finder<array_t, vector_t, tuple_t,
                                                  jagged_vector_t>;

    /// Surface finder definition
    // TODO: Move to volume builder
    using surfaces_finder_type =
        typename metadata::template surface_finder<array_t, vector_t, tuple_t,
                                                   jagged_vector_t>;

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
          _volume_finder(
              std::move(typename volume_finder::axis_p0_type{resource}),
              std::move(typename volume_finder::axis_p1_type{resource}),
              resource),
          _surfaces_finder(resource),
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
          _volume_finder(det_data._volume_finder_view),
          _surfaces_finder(det_data._surfaces_finder_view) {}

    /** Add a new volume and retrieve a reference to it
     *
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    DETRAY_HOST
    volume_type &new_volume(const array_t<scalar, 6> &bounds,
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
        dindex volume_index = _volume_finder.bin(p2);
        return _volumes[volume_index];
    }

    /** @return all surfaces - const access */
    DETRAY_HOST_DEVICE
    inline const auto &surfaces() const { return _surfaces; }

    /** @return all surfaces - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &surfaces() { return _surfaces; }

    /** @return a surface by index - const access */
    DETRAY_HOST_DEVICE
    inline const auto &surface_by_index(dindex sfidx) const {
        return _surfaces[sfidx];
    }

    /** @return a surface by index - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &surface_by_index(dindex sfidx) { return _surfaces[sfidx]; }

    /** @return all surface/portal masks in the geometry - const access */
    DETRAY_HOST_DEVICE
    inline auto &mask_store() const { return _masks; }

    /** @return all surface/portal masks in the geometry - non-const access */
    DETRAY_HOST_DEVICE
    inline auto &mask_store() { return _masks; }

    /** Add pre-built mask store
     *
     * @param masks the conatiner for surface masks
     */
    DETRAY_HOST
    inline void add_mask_store(mask_container &&msks) {
        _masks = std::move(msks);
    }

    /** Get all transform in an index range from the detector
     *
     * @param range The range of surfaces/portals in the transform store
     * @param ctx The context of the call
     *
     * @return ranged iterator to the object transforms
     */
    DETRAY_HOST_DEVICE
    inline auto transform_store(const dindex_range &range,
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
    inline const auto &transform_store(const context & /*ctx*/ = {}) const {
        return _transforms;
    }

    DETRAY_HOST_DEVICE
    inline auto &transform_store(const context & /*ctx*/ = {}) {
        return _transforms;
    }

    /** Add pre-built transform store
     *
     * @param transf the constianer for surface transforms
     */
    DETRAY_HOST
    inline void add_transform_store(transform_container &&transf) {
        _transforms = std::move(transf);
    }

    /** Get all available data from the detector without std::tie
     *
     * @param ctx The context of the call
     *
     * @return a struct that contains references to all relevant containers.
     */
    DETRAY_HOST_DEVICE
    auto data(const context & /*ctx*/ = {}) const {
        struct data_core {
            const dvector<volume_type> &volumes;
            const transform_container &transforms;
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
        detector_components &&...components) noexcept(false) {
        // Fill according to type, starting at type '0' (see 'masks')
        fill_containers(ctx, std::forward<detector_components>(components)...);
    }

    template <typename grid_type>
    DETRAY_HOST inline void add_surfaces_grid(const context ctx,
                                              volume_type &vol,
                                              grid_type &surfaces_grid) {
        // iterate over surfaces to fill the grid
        for (const auto [surf_idx, surf] : enumerate(_surfaces, vol)) {
            if (surf.get_grid_status() == true) {
                auto sidx = surf_idx;

                auto &trf =
                    _transforms.contextual_transform(ctx, surf.transform());
                auto tsl = trf.translation();

                if (vol.get_grid_type() ==
                    volume_type::grid_type::e_z_phi_grid) {

                    point2 location{tsl[2], algebra::getter::phi(tsl)};
                    surfaces_grid.populate(location, std::move(sidx));

                } else if (vol.get_grid_type() ==
                           volume_type::grid_type::e_r_phi_grid) {

                    point2 location{algebra::getter::perp(tsl),
                                    algebra::getter::phi(tsl)};
                    surfaces_grid.populate(location, std::move(sidx));
                }
            }
        }

        // add surfaces grid into surfaces finder
        auto n_grids = _surfaces_finder.effective_size();
        _surfaces_finder[n_grids] = surfaces_grid;
        vol.set_surfaces_finder(n_grids);
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
        mask_container &msks,
        transform_filling_container &trfs) noexcept(false) {

        // Get the surfaces/portals for a mask type
        auto &typed_surfaces = surfaces[current_type];
        // Get the corresponding transforms
        const auto &object_transforms = trfs[current_type];
        // and the corresponding masks
        auto &object_masks =
            msks.template group<mask_container::to_id(current_type)>();

        if (not object_transforms.empty(ctx) and not typed_surfaces.empty()) {
            // Current offsets into detectors containers
            const auto trsf_offset = _transforms.size(ctx);
            const auto mask_offset =
                _masks.template size<mask_container::to_id(current_type)>();

            // Fill the correct mask type
            _masks.add_masks(object_masks);
            _transforms.append(ctx, std::move(std::get<current_type>(trfs)));

            // Update the surfaces mask link
            for (auto &obj : typed_surfaces) {
                obj.update_mask(mask_offset);
                obj.update_transform(trsf_offset);
            }

            // Now put the updated objects into the geometry
            const auto offset = _surfaces.size();
            _surfaces.reserve(_surfaces.size() + typed_surfaces.size());
            _surfaces.insert(_surfaces.end(), typed_surfaces.begin(),
                             typed_surfaces.end());

            volume.update_range({offset, _surfaces.size()});
        }

        // Next mask type
        if constexpr (current_type <
                      std::tuple_size_v<typename mask_container::mask_tuple> -
                          1) {
            return fill_containers<current_type + 1, surface_container>(
                ctx, volume, surfaces, msks, trfs);
        }
        // update n_max_objects_per_volume
        else {
            _n_max_objects_per_volume =
                std::max(_n_max_objects_per_volume, volume.n_objects());
        }

        // If no mask type fits, don't fill the data.
    }

    /** Add the volume grid - move semantics
     *
     * @param v_grid the volume grid to be added
     */
    DETRAY_HOST
    inline void add_volume_finder(volume_finder &&v_grid) {
        _volume_finder = std::move(v_grid);
    }

    /** @return the volume grid - const access */
    DETRAY_HOST_DEVICE
    inline const volume_finder &volume_search_grid() const {
        return _volume_finder;
    }

    DETRAY_HOST_DEVICE
    inline volume_finder &volume_search_grid() { return _volume_finder; }

    DETRAY_HOST_DEVICE
    inline const surfaces_finder_type &get_surfaces_finder() const {
        return _surfaces_finder;
    }

    DETRAY_HOST_DEVICE
    inline surfaces_finder_type &get_surfaces_finder() {
        return _surfaces_finder;
    }

    DETRAY_HOST_DEVICE
    inline dindex get_n_max_objects_per_volume() const {
        return _n_max_objects_per_volume;
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
               << v.template n_objects<objects::e_surface>() << " surfaces "
               << std::endl;

            ss << "                 "
               << v.template n_objects<objects::e_portal>() << " portals "
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
    auto resource() const { return _resource; }

    private:
    /** Contains the geometrical relations */
    vector_t<volume_type> _volumes;

    /** All surfaces and portals in the geometry in contiguous memory */
    surface_container _surfaces;

    /** Keeps all of the transform data in contiguous memory*/
    transform_container _transforms;

    /** Surface and portal masks of the detector in contiguous memory */
    mask_container _masks;

    volume_finder _volume_finder;

    /* TODO: surfaces_finder needs to be refactored */
    surfaces_finder_type _surfaces_finder;

    vecmem::memory_resource *_resource = nullptr;

    // maximum number of surfaces per volume for navigation kernel candidates
    dindex _n_max_objects_per_volume = 0;
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
    using transform_container_t = typename detector_type::transform_container;
    using volume_finder_t = typename detector_type::volume_finder;
    using surfaces_finder_t = typename detector_type::surfaces_finder_type;

    detector_data(detector_type &det)
        : _volumes_data(vecmem::get_data(det.volumes())),
          _surfaces_data(vecmem::get_data(det.surfaces())),
          _masks_data(get_data(det.mask_store())),
          _transforms_data(get_data(det.transform_store())),
          _volume_finder_data(
              get_data(det.volume_search_grid(), *det.resource())),
          _surfaces_finder_data(
              get_data(det.get_surfaces_finder(), *det.resource())) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    mask_store_data<mask_container_t> _masks_data;
    static_transform_store_data<transform_container_t> _transforms_data;
    grid2_data<volume_finder_t> _volume_finder_data;
    surfaces_finder_data<surfaces_finder_t> _surfaces_finder_data;
};

/** A static inplementation of detector view for device
 */
template <typename detector_type>
struct detector_view {

    // type definitions
    using volume_t = typename detector_type::volume_type;
    using surface_t = typename detector_type::surface_type;
    using mask_container_t = typename detector_type::mask_container;
    using transform_container_t = typename detector_type::transform_container;
    using volume_finder_t = typename detector_type::volume_finder;
    using surfaces_finder_t = typename detector_type::surfaces_finder_type;

    detector_view(detector_data<detector_type> &det_data)
        : _volumes_data(det_data._volumes_data),
          _surfaces_data(det_data._surfaces_data),
          _masks_data(det_data._masks_data),
          _transforms_data(det_data._transforms_data),
          _volume_finder_view(det_data._volume_finder_data),
          _surfaces_finder_view(det_data._surfaces_finder_data) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    mask_store_data<mask_container_t> _masks_data;
    static_transform_store_data<transform_container_t> _transforms_data;
    grid2_view<volume_finder_t> _volume_finder_view;
    surfaces_finder_view<surfaces_finder_t> _surfaces_finder_view;
};

/** stand alone function for detector_data get function
 **/
template <
    typename detector_registry, template <typename, std::size_t> class array_t,
    template <typename...> class tuple_t, template <typename...> class vector_t,
    template <typename...> class jagged_vector_t, typename source_link>
inline detector_data<detector<detector_registry, array_t, tuple_t, vector_t,
                              jagged_vector_t, source_link> >
get_data(detector<detector_registry, array_t, tuple_t, vector_t,
                  jagged_vector_t, source_link> &det) {
    return det;
}

}  // namespace detray
