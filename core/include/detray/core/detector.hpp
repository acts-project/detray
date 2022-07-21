/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/tuple_vector_container.hpp"
#include "detray/core/detector_kernel.hpp"
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
#include "detray/intersection/intersection.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace detray {

/// @brief The detector definition.
///
/// This class is a heavily templated container aggregation, that owns all data
/// and sets the interface between geometry, navigator and surface finder
/// structures. Its view type is used to move the data between host and device.
///
/// @tparam metadata helper that defines collection and link types centrally
/// @tparam array_t the type of the internal array, must have STL semantics
/// @tparam tuple_t the type of the internal tuple, must have STL semantics
/// @tparam vector_t the type of the internal array, must have STL semantics
/// @tparam source_link the surface source link
template <typename metadata,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          typename source_link = dindex>
class detector {

    public:
    /// Algebra types
    using scalar_type = scalar;

    using point3 = __plugin::point3<scalar_type>;
    using vector3 = __plugin::vector3<scalar_type>;
    using point2 = __plugin::point2<scalar_type>;

    // Raw container types
    template <typename T>
    using vector_type = vector_t<T>;

    /// In case the detector needs to be printed
    using name_map = std::map<dindex, std::string>;

    /// Forward the alignable transform container (surface placements) and
    /// the geo context (e.g. for alignment)
    using transform_container =
        typename metadata::template transform_store<vector_t>;
    using transform3 = typename transform_container::transform3;
    using transform_link = typename transform_container::link_type;
    using context = typename transform_container::context;

    /// Forward mask types that are present in this detector
    using masks = typename metadata::mask_definitions;
    using mask_container =
        typename masks::template store_type<tuple_t, vector_t>;

    /// Forward material types that are present in this detector
    using materials = typename metadata::material_definitions;
    using material_container =
        typename materials::template store_type<tuple_t, vector_t>;

    /// Surface Finders: structures that enable neigborhood searches in the
    /// detector geometry during navigation. Can be different in each volume
    using sf_finders = typename metadata::template sf_finder_definitions<
        array_t, vector_t, tuple_t, jagged_vector_t>;
    using sf_finder_container =
        typename sf_finders::template store_type<tuple_t, array_t>;

    // TODO: Move to the following to volume builder

    /// The surface takes a mask (defines the local coordinates and the surface
    /// extent), its material, a link to an element in the transform container
    /// to define its placement and a source link to the object it represents.
    using surface_type = surface<masks, materials, transform_link, source_link>;
    /// Define the different kinds of surfaces that are present in the detector
    /// Can model the distinction between portals and sensitive surfaces
    using objects =
        typename metadata::template object_definitions<surface_type>;
    using surface_container = vector_t<surface_type>;
    /// Volume type
    using volume_type =
        volume<objects, scalar_type, typename sf_finders::link_type, array_t>;

    /// Volume finder definition: Make volume index available from track
    /// position
    using volume_finder =
        typename metadata::template volume_finder<array_t, vector_t, tuple_t,
                                                  jagged_vector_t>;

    /// Surface finder definition: Make neighboring surfaces available from
    /// track position.
    using surfaces_finder_type =
        typename metadata::template surface_finder<array_t, vector_t, tuple_t,
                                                   jagged_vector_t>;

    detector() = delete;

    /// Allowed costructor
    /// @param resource memory resource for the allocation of members
    DETRAY_HOST
    detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _sf_finders(resource),
          _volume_finder(
              std::move(typename volume_finder::axis_p0_type{resource}),
              std::move(typename volume_finder::axis_p1_type{resource}),
              resource),
          _surfaces_finder(resource),
          _resource(&resource) {}

    /// Constructor with detector_data
    template <typename detector_data_type,
              std::enable_if_t<!std::is_base_of_v<vecmem::memory_resource,
                                                  detector_data_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE detector(detector_data_type &det_data)
        : _volumes(det_data._volumes_data),
          _surfaces(det_data._surfaces_data),
          _transforms(det_data._transforms_data),
          _masks(det_data._masks_data),
          _materials(det_data._materials_data),
          _volume_finder(det_data._volume_finder_view),
          _surfaces_finder(det_data._surfaces_finder_view) {}

    /// Add a new volume and retrieve a reference to it
    ///
    /// @param bounds of the volume, they are expected to be already attaching
    /// @param surfaces_finder_entry of the volume, where to entry the surface
    /// finder
    ///
    /// @return non-const reference to the new volume
    DETRAY_HOST
    volume_type &new_volume(
        const array_t<scalar, 6> &bounds,
        typename volume_type::sf_finder_link_type sf_finder_link = {
            sf_finders::id::e_brute_force, dindex_invalid}) {
        volume_type &cvolume = _volumes.emplace_back(bounds);
        cvolume.set_index(_volumes.size() - 1);
        cvolume.set_sf_finder(sf_finder_link);

        return cvolume;
    }

    /// @return the sub-volumes of the detector - const access
    DETRAY_HOST_DEVICE
    inline auto volumes() const -> const vector_t<volume_type> & {
        return _volumes;
    }

    /// @return the sub-volumes of the detector - non-const access
    DETRAY_HOST_DEVICE
    inline auto volumes() -> vector_t<volume_type> & { return _volumes; }

    /// @return the volume by @param volume_index - const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_index(dindex volume_index) const
        -> const volume_type & {
        return _volumes[volume_index];
    }

    /// @return the volume by @param volume_index - non-const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_index(dindex volume_index) -> volume_type & {
        return _volumes[volume_index];
    }

    /// @return the volume by @param position - const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_pos(const point3 &p) const -> const volume_type & {
        point2 p2 = {getter::perp(p), p[2]};
        dindex volume_index = _volume_finder.bin(p2);
        return _volumes[volume_index];
    }

    /// @return all surfaces - const access
    DETRAY_HOST_DEVICE
    inline auto surfaces() const -> const surface_container & {
        return _surfaces;
    }

    /// @return all surfaces - non-const access
    /*DETRAY_HOST_DEVICE
    inline auto &surfaces() { return _surfaces; }*/

    /// @return a surface by index - const access
    DETRAY_HOST_DEVICE
    inline auto surface_by_index(dindex sfidx) const -> const surface_type & {
        return _surfaces[sfidx];
    }

    /// @return a surface by index - non-const access
    /*DETRAY_HOST_DEVICE
    inline auto &surface_by_index(dindex sfidx) { return _surfaces[sfidx]; }*/

    /// @return all surface/portal masks in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto mask_store() const -> const mask_container & { return _masks; }

    /// @return all surface/portal masks in the geometry - non-const access
    /*DETRAY_HOST_DEVICE
    inline auto &mask_store() { return _masks; }*/

    /// @return all materials in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto material_store() const -> const material_container & {
        return _materials;
    }

    /// @return all materials in the geometry - non-const access
    /*DETRAY_HOST_DEVICE
    inline auto &material_store() { return _materials; }*/

    /// Add pre-built mask store
    ///
    /// @param masks the conatiner for surface masks
    DETRAY_HOST
    inline void add_mask_store(mask_container &&msks) {
        _masks = std::move(msks);
    }

    /// Get all transform in an index range from the detector
    ///
    /// @param range The range of surfaces in the transform store
    /// @param ctx The context of the call
    ///
    /// @return ranged iterator to the surface transforms
    DETRAY_HOST_DEVICE
    inline auto transform_store(const dindex_range &range,
                                const context &ctx = {}) const {
        return _transforms.range(range, ctx);
    }

    /// Get all transform in an index range from the detector - const
    ///
    /// @param ctx The context of the call
    ///
    /// @return detector transform store
    DETRAY_HOST_DEVICE
    inline auto transform_store(const context & /*ctx*/ = {}) const
        -> const transform_container & {
        return _transforms;
    }

    /*DETRAY_HOST_DEVICE
    inline auto &transform_store(const context & /*ctx*/ //= {}) {
    //    return _transforms;
    //}*/

    /// Add pre-built transform store
    ///
    /// @param transf the constianer for surface transforms
    DETRAY_HOST
    inline void add_transform_store(transform_container &&transf) {
        _transforms = std::move(transf);
    }

    /// Get all available data from the detector without std::tie
    ///
    /// @param ctx The context of the call
    ///
    /// @return a struct that contains references to all relevant containers.
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

    // TODO: Provide grid builder structure separate from the detector
    template <typename sf_finder_t>
    DETRAY_HOST void add_sf_finder(const context ctx, volume_type &vol,
                                   sf_finder_t &surface_finder) {
        // Get id for the surface finder link in the volume
        constexpr typename sf_finders::id sf_finder_id =
            sf_finders::template get_id<sf_finder_t>();

        // iterate over surfaces to fill the grid
        for (const auto [surf_idx, surf] : enumerate(_surfaces, vol)) {
            if (surf.get_grid_status() == true) {
                dindex sidx = surf_idx;

                const transform3 &trf =
                    _transforms.contextual_transform(ctx, surf.transform());
                auto tsl = trf.translation();

                if constexpr (sf_finder_id == sf_finders::id::e_z_phi_grid) {
                    std::cout << "Fill z_phi grid for vol " << vol.index()
                              << std::endl;

                    point2 location{tsl[2], algebra::getter::phi(tsl)};
                    surface_finder.populate(location, std::move(sidx));

                } else if constexpr (sf_finder_id ==
                                     sf_finders::id::e_r_phi_grid) {
                    std::cout << "Fill r_phi grid for vol " << vol.index()
                              << std::endl;

                    point2 location{algebra::getter::perp(tsl),
                                    algebra::getter::phi(tsl)};
                    surface_finder.populate(location, std::move(sidx));
                }
            }
        }

        // add surfaces grid into surfaces finder container
        auto &sf_finder_group = _sf_finders.template group<sf_finder_id>();
        // Find index of this surface finder
        std::size_t sf_finder_idx = 0;
        for (unsigned int i_s = 0; i_s < sf_finder_group.size(); i_s++) {
            if (!sf_finder_group[i_s].data().empty()) {
                sf_finder_idx++;
            }
        }
        sf_finder_group.at(sf_finder_idx) = surface_finder;
        vol.set_sf_finder(sf_finder_id, sf_finder_idx);
    }

    /// Add a new full set of detector components (e.g. transforms or volumes)
    /// according to given context.
    ///
    /// @param ctx is the context of the call
    /// @param vol is the target volume
    /// @param surfaces_per_vol is the surface vector per volume
    /// @param masks_per_vol is the mask container per volume
    /// @param materials_per_vol is the material container per volume
    /// @param trfs_per_vol is the transform vector per volume
    ///
    /// @note can throw an exception if input data is inconsistent
    // TODO: Provide volume builder structure separate from the detector
    DETRAY_HOST void add_objects_per_volume(
        const context ctx, volume_type &vol,
        surface_container &surfaces_per_vol, mask_container &masks_per_vol,
        material_container &materials_per_vol,
        transform_container &trfs_per_vol) noexcept(false) {

        // Append transforms
        const auto trf_offset = _transforms.size(ctx);
        _transforms.append(ctx, std::move(trfs_per_vol));

        // Update mask, material and transform index of surfaces
        for (auto &sf : surfaces_per_vol) {
            _masks.template execute<mask_index_update>(sf.mask_type(), sf);
            _materials.template execute<material_index_update>(
                sf.material_type(), sf);
            sf.update_transform(trf_offset);
        }

        // Append surfaces
        const auto sf_offset = _surfaces.size();
        _surfaces.reserve(sf_offset + surfaces_per_vol.size());
        _surfaces.insert(_surfaces.end(), surfaces_per_vol.begin(),
                         surfaces_per_vol.end());

        // Update the surface range per volume
        vol.update_range({sf_offset, _surfaces.size()});

        // Append mask and material container
        _masks.append_container(masks_per_vol);
        _materials.append_container(materials_per_vol);

        // Update max objects per volume
        _n_max_objects_per_volume =
            std::max(_n_max_objects_per_volume, vol.n_objects());
    }

    /// Add the volume grid - move semantics
    ///
    /// @param v_grid the volume grid to be added
    DETRAY_HOST
    inline void add_volume_finder(volume_finder &&v_grid) {
        _volume_finder = std::move(v_grid);
    }

    /// @return the volume grid - const access
    DETRAY_HOST_DEVICE
    inline const volume_finder &volume_search_grid() const {
        return _volume_finder;
    }

    /// @returns const access to the detector's volume search structure
    /*DETRAY_HOST_DEVICE
    inline volume_finder &volume_search_grid() { return _volume_finder; }*/

    /// @returns access to the surface finder container - non-const access
    // TODO: remove once possible
    DETRAY_HOST_DEVICE
    inline auto sf_finder() -> surfaces_finder_type & {
        return _surfaces_finder;
    }

    /// @returns access to the surface finder container - non-const access
    // TODO: remove once possible
    DETRAY_HOST_DEVICE
    inline auto sf_finder_store() -> sf_finder_container & {
        return _sf_finders;
    }

    /// @returns access to the surface finder container
    DETRAY_HOST_DEVICE
    inline auto sf_finder_store() const -> const sf_finder_container & {
        return _sf_finders;
    }

    /// @returns the maximum number of surfaces (sensitive + portal) in all
    /// volumes.
    DETRAY_HOST_DEVICE
    inline dindex get_n_max_objects_per_volume() const {
        return _n_max_objects_per_volume;
    }

    /// @param names maps a volume to its string representation.
    /// @returns a string representation of the detector.
    DETRAY_HOST
    const std::string to_string(const name_map &names) const {
        std::stringstream ss;

        ss << "[>] Detector '" << names.at(0) << "' has " << _volumes.size()
           << " volumes." << std::endl;
        ss << "    contains  " << _surfaces_finder.size()
           << " local surface finders." << std::endl;

        for (const auto &[i, v] : enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << ": " << std::endl;
            ss << " - name: '" << v.name(names) << "'" << std::endl;

            ss << "     contains    "
               << v.template n_objects<objects::e_surface>() << " surfaces "
               << std::endl;

            ss << "                 "
               << v.template n_objects<objects::e_portal>() << " portals "
               << std::endl;

            if (v.sf_finder_index() != dindex_invalid) {
                ss << "  sf finder id " << v.sf_finder_type()
                   << "  sf finders idx " << v.sf_finder_index() << std::endl;
            }

            const auto &bounds = v.bounds();
            ss << "     bounds r = (" << bounds[0] << ", " << bounds[1] << ")"
               << std::endl;
            ss << "            z = (" << bounds[2] << ", " << bounds[3] << ")"
               << std::endl;
        }

        return ss.str();
    }

    /// @return the pointer of memoery resource - non-const access
    DETRAY_HOST
    auto resource() -> vecmem::memory_resource * { return _resource; }

    /// @return the pointer of memoery resource
    DETRAY_HOST
    auto resource() const -> const vecmem::memory_resource * {
        return _resource;
    }

    private:
    /// Contains the detector sub-volumes.
    vector_t<volume_type> _volumes;

    /// All surfaces (sensitive and portal) in the geometry in contiguous memory
    surface_container _surfaces;

    /// Keeps all of the transform data in contiguous memory
    transform_container _transforms;

    /// Masks of all surfaces in the geometry in contiguous memory
    mask_container _masks;

    /// Materials of all surfaces in the geometry in contiguous memory
    material_container _materials;

    /// All surface finder data structures that are used in the detector volumes
    sf_finder_container _sf_finders;

    /// Search structure for volumes
    volume_finder _volume_finder;

    /// Container for all surface search structures per volume, e.g. grids
    // TODO: surfaces_finder needs to be refactored
    surfaces_finder_type _surfaces_finder;

    /// The memory resource represents how and where (host, device, managed)
    /// the memory for the detector containers is allocated
    vecmem::memory_resource *_resource = nullptr;

    /// Maximum number of surfaces per volume. Used to estimate size of
    /// candidates vector in the navigator. Is determined during detector
    /// building.
    dindex _n_max_objects_per_volume = 0;
};

/// @brief A static inplementation of detector data for device
template <typename detector_type>
struct detector_data {

    // type definitions
    using volume_t = typename detector_type::volume_type;
    using surface_t = typename detector_type::surface_type;
    using mask_container_t = typename detector_type::mask_container;
    using material_container_t = typename detector_type::material_container;
    using transform_container_t = typename detector_type::transform_container;
    using volume_finder_t = typename detector_type::volume_finder;
    using surfaces_finder_t = typename detector_type::surfaces_finder_type;

    detector_data(detector_type &det)
        : _volumes_data(vecmem::get_data(det.volumes())),
          _surfaces_data(vecmem::get_data(det.surfaces())),
          _masks_data(get_data(det.mask_store())),
          _materials_data(get_data(det.material_store())),
          _transforms_data(get_data(det.transform_store())),
          _volume_finder_data(
              get_data(det.volume_search_grid(), *det.resource())),
          _surfaces_finder_data(
              get_data(det.get_surfaces_finder(), *det.resource())) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    tuple_vector_container_data<mask_container_t> _masks_data;
    tuple_vector_container_data<material_container_t> _materials_data;
    static_transform_store_data<transform_container_t> _transforms_data;
    grid2_data<volume_finder_t> _volume_finder_data;
    surfaces_finder_data<surfaces_finder_t> _surfaces_finder_data;
};

/// @brief A static inplementation of detector view for device
template <typename detector_type>
struct detector_view {

    // type definitions
    using volume_t = typename detector_type::volume_type;
    using surface_t = typename detector_type::surface_type;
    using mask_container_t = typename detector_type::mask_container;
    using material_container_t = typename detector_type::material_container;
    using transform_container_t = typename detector_type::transform_container;
    using volume_finder_t = typename detector_type::volume_finder;
    using surfaces_finder_t = typename detector_type::surfaces_finder_type;

    detector_view(detector_data<detector_type> &det_data)
        : _volumes_data(det_data._volumes_data),
          _surfaces_data(det_data._surfaces_data),
          _masks_data(det_data._masks_data),
          _materials_data(det_data._materials_data),
          _transforms_data(det_data._transforms_data),
          _volume_finder_view(det_data._volume_finder_data),
          _surfaces_finder_view(det_data._surfaces_finder_data) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    tuple_vector_container_data<mask_container_t> _masks_data;
    tuple_vector_container_data<material_container_t> _materials_data;
    static_transform_store_data<transform_container_t> _transforms_data;
    grid2_view<volume_finder_t> _volume_finder_view;
    surfaces_finder_view<surfaces_finder_t> _surfaces_finder_view;
};

/// stand alone function for that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename metadata, template <typename, std::size_t> class array_t,
          template <typename...> class tuple_t,
          template <typename...> class vector_t,
          template <typename...> class jagged_vector_t, typename source_link>
inline detector_data<detector<metadata, array_t, tuple_t, vector_t,
                              jagged_vector_t, source_link> >
get_data(detector<metadata, array_t, tuple_t, vector_t, jagged_vector_t,
                  source_link> &det) {
    return det;
}

}  // namespace detray
