/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/detector_kernel.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>

// System include(s)
#include <map>
#include <sstream>
#include <string>
#include <iostream>

namespace detray {
    using _toy_bfield = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::float3>>>>;


/// @brief Forward declaration of a detector view type
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
struct detector_view;

/// @brief The detector definition.
///
/// This class is a heavily templated container aggregation, that owns all data
/// and sets the interface between geometry, navigator and surface finder
/// structures. Its view type is used to move the data between host and device.
///
/// @tparam metadata helper that defines collection and link types centrally
/// @tparam bfield_t the type of the b-field frontend
/// @tparam container_t type collection of the underlying containers
/// @tparam source_link the surface source link
template <typename metadata, template <typename> class bfield_t = covfie::field,
          typename container_t = host_container_types>
class detector {

    // Allow the building of the detector containers
    friend class volume_builder_interface<
        detector<metadata, bfield_t, container_t>>;

    public:
    /// Algebra types
    /// @TODO: scalar as a template parameter
    using scalar_type = scalar;

    using point3 = __plugin::point3<scalar_type>;
    using vector3 = __plugin::vector3<scalar_type>;
    using point2 = __plugin::point2<scalar_type>;

    using bfield_backend_type = typename metadata::bfield_backend_t;

    using bfield_type = bfield_t<bfield_backend_type>;

    /// Raw container types
    template <typename T, std::size_t N>
    using array_type = typename container_t::template array_type<T, N>;
    template <typename T>
    using vector_type = typename container_t::template vector_type<T>;
    template <typename... T>
    using tuple_type = typename container_t::template tuple_type<T...>;
    template <typename T>
    using jagged_vector_type =
        typename container_t::template jagged_vector_type<T>;

    /// In case the detector needs to be printed
    using name_map = std::map<dindex, std::string>;

    /// Forward the alignable transform container (surface placements) and
    /// the geo context (e.g. for alignment)
    using transform_container =
        typename metadata::template transform_store<vector_type>;
    using transform3 = typename transform_container::value_type;
    using transform_link = typename transform_container::link_type;
    using geometry_context = typename transform_container::context_type;

    /// Forward mask types that are present in this detector
    using mask_container =
        typename metadata::template mask_store<tuple_type, vector_type>;
    using masks = typename mask_container::value_types;
    using mask_link = typename mask_container::single_link;

    /// Forward mask types that are present in this detector
    using material_container =
        typename metadata::template material_store<tuple_type, vector_type>;
    using materials = typename material_container::value_types;
    using material_link = typename material_container::single_link;

    /// Surface Finders: structures that enable neigborhood searches in the
    /// detector geometry during navigation. Can be different in each volume
    using sf_finder_container =
        typename metadata::template surface_finder_store<tuple_type,
                                                         container_t>;
    using sf_finders = typename sf_finder_container::value_types;
    using sf_finder_link = typename sf_finder_container::single_link;

    // TODO: Move to the following to volume builder

    /// The surface takes a mask (defines the local coordinates and the surface
    /// extent), its material, a link to an element in the transform container
    /// to define its placement and a source link to the object it represents.
    using surface_type = typename metadata::surface_type;

    using surface_container = vector_type<surface_type>;
    /// Volume type
    using geo_obj_ids = typename metadata::geo_objects;
    using volume_type =
        detector_volume<geo_obj_ids, typename metadata::object_link_type,
                        sf_finder_link, scalar_type>;

    /// Volume finder definition: Make volume index available from track
    /// position
    using volume_finder =
        typename metadata::template volume_finder<container_t>;

    using detector_view_type =
        detector_view<metadata, covfie::field, host_container_types>;

    detector() = delete;

    /// Allowed costructor
    /// @param resource memory resource for the allocation of members
    DETRAY_HOST
    detector(vecmem::memory_resource &resource, bfield_type &&field)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _sf_finders(resource),
          _volume_finder(resource),
          _resource(&resource),
          _bfield(field) {}

    static covfie::field<_toy_bfield> load_toy_field() {
        std::ifstream stream(std::getenv("DETRAY_BFIELD_FILE"), std::ifstream::binary);

        if (!stream.good()) {
            std::cout << "Lol the file failed to load!!" << std::endl;
        }

        covfie::field<_toy_bfield> field(stream);

        stream.close();

        return field;
    }

    /// Constructor with simplified constant-zero B-field
    /// @param resource memory resource for the allocation of members
    template<typename toy_bfield = _toy_bfield, std::enable_if_t<std::is_same_v<typename metadata::bfield_backend_t, toy_bfield>, bool> = true>
    DETRAY_HOST
    detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _sf_finders(resource),
          _volume_finder(resource),
          _resource(&resource),
          _bfield(load_toy_field()) {}

    template<typename toy_bfield = _toy_bfield, std::enable_if_t<!std::is_same_v<typename metadata::bfield_backend_t, toy_bfield>, bool> = true>
    DETRAY_HOST
    detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _sf_finders(resource),
          _volume_finder(resource),
          _resource(&resource),
          _bfield(typename bfield_type::backend_t::configuration_t{0.f, 0.f,
                                                                   0.f}) {}

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
          _sf_finders(det_data._sf_finder_data),
          _volume_finder(det_data._volume_finder_data),
          _bfield(det_data._bfield_view) {}

    /// Add a new volume and retrieve a reference to it
    ///
    /// @param bounds of the volume, they are expected to be already attaching
    /// @param sf_finder_link of the volume, where to entry the surface finder
    ///
    /// @return non-const reference to the new volume
    DETRAY_HOST
    volume_type &new_volume(
        const volume_id id, const array_type<scalar, 6> &bounds,
        typename volume_type::sf_finder_link_type srf_finder_link = {
            sf_finders::id::e_default, dindex_invalid}) {
        volume_type &cvolume = _volumes.emplace_back(id, bounds);
        cvolume.set_index(_volumes.size() - 1);
        cvolume.set_sf_finder(srf_finder_link);

        return cvolume;
    }

    /// @return the sub-volumes of the detector - const access
    DETRAY_HOST_DEVICE
    inline auto volumes() const -> const vector_type<volume_type> & {
        return _volumes;
    }

    /// @return the sub-volumes of the detector - non-const access
    DETRAY_HOST_DEVICE
    inline auto volumes() -> vector_type<volume_type> & { return _volumes; }

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
        // Grid with replace-populator: only one entry per bin
        dindex volume_index = *_volume_finder.search(p);
        return _volumes[volume_index];
    }

    /// @return all surfaces - const access
    DETRAY_HOST_DEVICE
    inline auto surfaces() const -> const surface_container & {
        return _surfaces;
    }

    /// @return all surfaces - non-const access
    DETRAY_HOST_DEVICE
    inline auto surfaces() -> surface_container & { return _surfaces; }

    /// @return a surface by index - const access
    DETRAY_HOST_DEVICE
    inline auto surface_by_index(dindex sfidx) const -> const surface_type & {
        return _surfaces[sfidx];
    }

    /// @return a surface by index - non-const access
    DETRAY_HOST_DEVICE
    inline auto surface_by_index(dindex sfidx) -> surface_type & {
        return _surfaces[sfidx];
    }

    /// Append new surfaces to the detector
    DETRAY_HOST
    inline void append_surfaces(surface_container &&new_surfaces) {
        _surfaces.insert(_surfaces.end(), new_surfaces.begin(),
                         new_surfaces.end());
    }

    /// @return all surface/portal masks in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto mask_store() const -> const mask_container & { return _masks; }

    /// @return all surface/portal masks in the geometry - non-const access
    DETRAY_HOST_DEVICE
    inline auto mask_store() -> mask_container & { return _masks; }

    /// Append a new mask store to the detector
    DETRAY_HOST
    inline void append_masks(mask_container &&new_masks) {
        _masks.append(std::move(new_masks));
    }

    /// @return all materials in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto material_store() const -> const material_container & {
        return _materials;
    }

    /// @return all materials in the geometry - non-const access
    DETRAY_HOST_DEVICE
    inline auto material_store() -> material_container & { return _materials; }

    /// Append a new material store to the detector
    DETRAY_HOST
    inline void append_materials(material_container &&new_materials) {
        _materials.append(std::move(new_materials));
    }

    /// Get all transform in an index range from the detector - const
    ///
    /// @param ctx The context of the call
    ///
    /// @return detector transform store
    DETRAY_HOST_DEVICE
    inline auto transform_store(const geometry_context & /*ctx*/ = {}) const
        -> const transform_container & {
        return _transforms;
    }

    DETRAY_HOST_DEVICE
    inline auto transform_store(const geometry_context & /*ctx*/ = {})
        -> transform_container & {
        return _transforms;
    }

    /// Append a new transform store to the detector
    DETRAY_HOST
    inline void append_transforms(transform_container &&new_transforms,
                                  const geometry_context ctx = {}) {
        _transforms.append(std::move(new_transforms), ctx);
    }

    /// Get all available data from the detector without std::tie
    ///
    /// @param ctx The context of the call
    ///
    /// @return a struct that contains references to all relevant containers.
    DETRAY_HOST_DEVICE
    auto data(const geometry_context & /*ctx*/ = {}) const {
        struct data_core {
            const dvector<volume_type> &volumes;
            const transform_container &transforms;
            const mask_container &masks;
            const surface_container &surfaces;
        };
        return data_core{_volumes, _transforms, _masks, _surfaces};
    }

    /// Add a new full set of detector components (e.g. transforms or volumes)
    /// according to given geometry_context.
    ///
    /// @param ctx is the geometry_context of the call
    /// @param vol is the target volume
    /// @param surfaces_per_vol is the surface vector per volume
    ///                         (either portals, sensitives or passives)
    /// @param masks_per_vol is the mask container per volume
    /// @param trfs_per_vol is the transform vector per volume
    ///
    /// @note can throw an exception if input data is inconsistent
    // TODO: Provide volume builder structure separate from the detector
    template <geo_obj_ids surface_id = geo_obj_ids::e_portal>
    DETRAY_HOST auto add_objects_per_volume(
        const geometry_context ctx, volume_type &vol,
        surface_container &surfaces_per_vol, mask_container &masks_per_vol,
        transform_container &trfs_per_vol) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = _transforms.size(ctx);
        _transforms.append(std::move(trfs_per_vol), ctx);

        // Update mask, material and transform index of surfaces
        for (auto &sf : surfaces_per_vol) {
            _masks.template call<detail::mask_index_update>(sf.mask(), sf);
            sf.update_transform(trf_offset);
        }

        // Append surfaces
        const auto sf_offset = _surfaces.size();
        _surfaces.reserve(sf_offset + surfaces_per_vol.size());
        _surfaces.insert(_surfaces.end(), surfaces_per_vol.begin(),
                         surfaces_per_vol.end());

        // Update the surface range per volume
        vol.template update_obj_link<surface_id>(sf_offset,
                                                 surfaces_per_vol.size());

        // Append mask and material container
        _masks.append(std::move(masks_per_vol));

        // Update max objects per volume
        _n_max_objects_per_volume =
            std::max(_n_max_objects_per_volume, vol.n_objects());
    }

    /// Add a new full set of detector components (e.g. transforms or volumes)
    /// according to given geometry_context.
    ///
    /// @param ctx is the geometry_context of the call
    /// @param vol is the target volume
    /// @param surfaces_per_vol is the surface vector per volume
    /// @param masks_per_vol is the mask container per volume
    /// @param materials_per_vol is the material container per volume
    /// @param trfs_per_vol is the transform vector per volume
    ///
    /// @note can throw an exception if input data is inconsistent
    // TODO: Provide volume builder structure separate from the detector
    DETRAY_HOST
    auto add_objects_per_volume(
        const geometry_context ctx, volume_type &vol,
        surface_container &surfaces_per_vol, mask_container &masks_per_vol,
        transform_container &trfs_per_vol,
        material_container &materials_per_vol) noexcept(false) -> void {

        // Update material index of surfaces
        for (auto &sf : surfaces_per_vol) {
            _materials.template call<detail::material_index_update>(
                sf.material(), sf);
        }
        _materials.append(std::move(materials_per_vol));

        add_objects_per_volume(ctx, vol, surfaces_per_vol, masks_per_vol,
                               trfs_per_vol);
    }

    /// Add the volume grid - move semantics
    ///
    /// @param v_grid the volume grid to be added
    DETRAY_HOST
    inline auto add_volume_finder(volume_finder &&v_grid) -> void {
        _volume_finder = std::move(v_grid);
    }

    /// @return the volume grid - const access
    DETRAY_HOST_DEVICE
    inline auto volume_search_grid() const -> const volume_finder & {
        return _volume_finder;
    }

    /// @returns const access to the detector's volume search structure
    DETRAY_HOST_DEVICE
    inline auto volume_search_grid() -> volume_finder & {
        return _volume_finder;
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

    /// @brief Updates the maximum number of surfaces (sensitive + passive +
    /// portal) over all volumes.
    DETRAY_HOST_DEVICE
    inline auto update_n_max_objects_per_volume(std::size_t n_surfaces)
        -> void {
        _n_max_objects_per_volume =
            std::max(_n_max_objects_per_volume, n_surfaces);
    }

    /// @returns the maximum number of surfaces (sensitive + passive + portal)
    /// over all volumes.
    DETRAY_HOST_DEVICE
    inline auto get_n_max_objects_per_volume() const -> dindex {
        return _n_max_objects_per_volume;
    }

    DETRAY_HOST_DEVICE
    inline const bfield_type &get_bfield() const { return _bfield; }

    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const dindex sf_idx, const point3 &pos,
                                  const vector3 &dir) const {
        const auto &sf = surface_by_index(sf_idx);
        const auto ret =
            _masks.template call<detail::global_to_local<transform3>>(
                sf.mask(), _transforms[sf.transform()], pos, dir);
        return ret;
    }

    /// @param names maps a volume to its string representation.
    /// @returns a string representation of the detector.
    DETRAY_HOST
    auto to_string(const name_map &names) const -> std::string {
        std::stringstream ss;

        ss << "[>] Detector '" << names.at(0) << "' has " << _volumes.size()
           << " volumes." << std::endl;
        ss << " local surface finders." << std::endl;

        for (const auto [i, v] : detray::views::enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << ": " << std::endl;
            ss << " - name: '" << v.name(names) << "'" << std::endl;

            ss << "     contains    "
               << v.template n_objects<geo_obj_ids::e_sensitive>()
               << " sensitive surfaces " << std::endl;

            ss << "                 "
               << v.template n_objects<geo_obj_ids::e_portal>() << " portals "
               << std::endl;

            ss << "                 " << _sf_finders.n_collections()
               << " surface finders " << std::endl;

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
    auto *resource() { return _resource; }

    /// @return the pointer of memoery resource
    DETRAY_HOST
    const auto *resource() const { return _resource; }

    private:
    /// Contains the detector sub-volumes.
    vector_type<volume_type> _volumes;

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

    /// The memory resource represents how and where (host, device, managed)
    /// the memory for the detector containers is allocated
    vecmem::memory_resource *_resource = nullptr;

    /// Maximum number of surfaces per volume. Used to estimate size of
    /// candidates vector in the navigator. Is determined during detector
    /// building.
    dindex _n_max_objects_per_volume = 0;

    /// Storage for magnetic field data
    bfield_type _bfield;
};

/// @brief A static inplementation of detector data for device
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
struct detector_view {

    using detector_type = detector<metadata, bfield_t, container_t>;

    detector_view(detector_type &det)
        : _volumes_data(vecmem::get_data(det.volumes())),
          _surfaces_data(vecmem::get_data(det.surfaces())),
          _masks_data(get_data(det.mask_store())),
          _materials_data(get_data(det.material_store())),
          _transforms_data(get_data(det.transform_store())),
          _sf_finder_data(get_data(det.sf_finder_store())),
          _volume_finder_data(get_data(det.volume_search_grid())),
          _bfield_view(det.get_bfield()) {}

    // members
    vecmem::data::vector_view<typename detector_type::volume_type>
        _volumes_data;
    vecmem::data::vector_view<typename detector_type::surface_type>
        _surfaces_data;
    typename detector_type::mask_container::view_type _masks_data;
    typename detector_type::material_container::view_type _materials_data;
    typename detector_type::transform_container::view_type _transforms_data;
    typename detector_type::sf_finder_container::view_type _sf_finder_data;
    typename detector_type::volume_finder::view_type _volume_finder_data;
    typename detector_type::bfield_type::view_t _bfield_view;
};

/// stand alone function for that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
inline detector_view<metadata, bfield_t, container_t> get_data(
    detector<metadata, bfield_t, container_t> &det) {
    return {det};
}

}  // namespace detray
