/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector_kernel.hpp"
#include "detray/core/surfaces_finder.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/bin_association.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
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
/// @tparam bfield_t the type of the b-field frontend
/// @tparam container_t type collection of the underlying containers
/// @tparam source_link the surface source link
template <typename metadata, template <typename> class bfield_t = covfie::field,
          typename container_t = host_container_types,
          typename source_link = dindex>
class detector {

    public:
    /// Algebra types
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
    using context = typename transform_container::context_type;

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
    /*using sf_finders = typename metadata::template sf_finder_definitions<
        array_t, vector_t, tuple_t, jagged_vector_t>;
    using sf_finder_container =
        typename sf_finders::template store_type<tuple_t, array_t>;*/
    using sf_finder_container =
        typename metadata::template surface_finder_store<tuple_type,
                                                         container_t>;
    using sf_finders = typename sf_finder_container::value_types;
    using sf_finder_link = typename sf_finder_container::single_link;

    // TODO: Move to the following to volume builder

    /// The surface takes a mask (defines the local coordinates and the surface
    /// extent), its material, a link to an element in the transform container
    /// to define its placement and a source link to the object it represents.
    using surface_type =
        surface<mask_link, material_link, transform_link, source_link>;

    using surface_container = vector_type<surface_type>;
    /// Volume type
    using geo_obj_ids = typename metadata::geo_objects;
    using volume_type =
        detector_volume<geo_obj_ids, typename metadata::object_link_type,
                        sf_finder_link, scalar_type>;

    /// Volume finder definition: Make volume index available from track
    /// position
    using volume_finder = typename metadata::template volume_finder<
        array_type, vector_type, tuple_type, jagged_vector_type>;

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
          _volume_finder(
              std::move(typename volume_finder::axis_p0_type{resource}),
              std::move(typename volume_finder::axis_p1_type{resource}),
              resource),
          _resource(&resource),
          _bfield(field) {}

    /// Constructor with simplified constant-zero B-field
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
          _volume_finder(det_data._volume_finder_view),
          _bfield(det_data._bfield_view) {}

    /// Add a new volume and retrieve a reference to it
    ///
    /// @param bounds of the volume, they are expected to be already attaching
    /// @param sf_finder_link of the volume, where to entry the surface finder
    ///
    /// @return non-const reference to the new volume
    DETRAY_HOST
    volume_type &new_volume(
        const array_type<scalar, 6> &bounds,
        typename volume_type::sf_finder_link_type srf_finder_link = {
            sf_finders::id::e_default, dindex_invalid}) {
        volume_type &cvolume = _volumes.emplace_back(bounds);
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

    /// @return all surface/portal masks in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto mask_store() const -> const mask_container & { return _masks; }

    /// @return all surface/portal masks in the geometry - non-const access
    DETRAY_HOST_DEVICE
    inline auto mask_store() -> mask_container & { return _masks; }

    /// Add pre-built mask store
    ///
    /// @param masks the conatiner for surface masks
    DETRAY_HOST
    inline void add_mask_store(mask_container &&msks) {
        _masks = std::move(msks);
    }

    /// @return all materials in the geometry - const access
    DETRAY_HOST_DEVICE
    inline auto material_store() const -> const material_container & {
        return _materials;
    }

    /// @return all materials in the geometry - non-const access
    DETRAY_HOST_DEVICE
    inline auto material_store() -> material_container & { return _materials; }

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

    DETRAY_HOST_DEVICE
    inline auto transform_store(const context & /*ctx*/ = {})
        -> transform_container & {
        return _transforms;
    }

    /// Add pre-built transform store
    ///
    /// @param transf the constianer for surface transforms
    DETRAY_HOST
    inline auto add_transform_store(transform_container &&transf) -> void {
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

    /// Add a grid to the surface finder store of the detector
    ///
    /// New surface finder id can be been given explicitly. That is helpful, if
    /// multiple sf finders have the same type in the tuple container. Otherwise
    /// it is determined automatically.
    ///
    /// @param vol the volume the surface finder should be added to
    /// @param grid the grid that should be added
    // TODO: Provide grid builder structure separate from the detector
    /*template <typename grid_t, typename sf_finders::id grid_id =
                                   sf_finders::template get_id<grid_t>()>
    DETRAY_HOST auto add_grid(volume_type &vol, grid_t &grid) -> void {

        // Add surfaces grid to surfaces finder container
        auto &grid_group = _sf_finders.template group<grid_id>();

        // Find correct index for this surface finder
        std::size_t grid_idx = 0;
        for (unsigned int i_s = 0; i_s < grid_group.size(); i_s++) {
            if (!grid_group.at(i_s).data().empty()) {
                grid_idx++;
            }
        }
        grid_group.at(grid_idx) = std::move(grid);
        vol.set_sf_finder(grid_id, grid_idx);
    }

    /// Fill a grid surface finder by bin association, then add it to the
    /// detector.
    ///
    /// New surface finder id can be been given explicitly. That is helpful, if
    /// multiple sf finders have the same type in the tuple container. Otherwise
    /// it is determined automatically.
    ///
    /// @param ctx the geometry context
    /// @param vol the volume the surface finder should be added to
    /// @param grid the grid that should be added
    // TODO: Provide grid builder structure separate from the detector
    template <typename grid_t>
    DETRAY_HOST auto fill_grid(const context ctx, volume_type &vol,
                               grid_t &grid) -> void {

        // Fill the volumes surfaces into the grid
        bin_association(ctx, *this, vol, grid, {0.1, 0.1}, false);
    }

    /// Detector interface to add surface finders
    ///
    /// @param ctx the geometry context
    /// @param vol the volume the surface finder should be added to
    /// @param grid the grid that should be added
    template <typename sf_finder_t,
              typename sf_finders::id sf_finder_id =
                  sf_finders::template get_id<sf_finder_t>()>
    DETRAY_HOST auto add_sf_finder(const context ctx, volume_type &vol,
                                   sf_finder_t &sf_finder) -> void {

        // For now, only implemented for grids
        fill_grid(ctx, vol, sf_finder);
        add_grid<sf_finder_t, sf_finder_id>(vol, sf_finder);
    }*/

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
    DETRAY_HOST
    auto add_objects_per_volume(
        const context ctx, volume_type &vol,
        surface_container &surfaces_per_vol, mask_container &masks_per_vol,
        material_container &materials_per_vol,
        transform_container &trfs_per_vol) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = _transforms.size(ctx);
        _transforms.append(std::move(trfs_per_vol), ctx);

        // Update mask, material and transform index of surfaces
        for (auto &sf : surfaces_per_vol) {
            _masks.template call<mask_index_update>(sf.mask(), sf);
            _materials.template call<material_index_update>(sf.material(), sf);
            sf.update_transform(trf_offset);
        }

        // Append surfaces
        const auto sf_offset = _surfaces.size();
        _surfaces.reserve(sf_offset + surfaces_per_vol.size());
        _surfaces.insert(_surfaces.end(), surfaces_per_vol.begin(),
                         surfaces_per_vol.end());

        // Update the surface range per volume
        vol.update_obj_link({sf_offset, _surfaces.size()});

        // Append mask and material container
        _masks.append(std::move(masks_per_vol));
        _materials.append(std::move(materials_per_vol));

        // Update max objects per volume
        _n_max_objects_per_volume =
            std::max(_n_max_objects_per_volume, vol.n_objects());
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

    /// @returns the maximum number of surfaces (sensitive + portal) in all
    /// volumes.
    DETRAY_HOST_DEVICE
    inline auto get_n_max_objects_per_volume() const -> dindex {
        return _n_max_objects_per_volume;
    }

    DETRAY_HOST_DEVICE
    inline const bfield_type &get_bfield() const { return _bfield; }

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

            ss << "                 " << /*_sf_finders.size(v.sf_finder_type())
               << " surface finders " <<*/
                std::endl;

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

    /** Storage for magnetic field data */
    bfield_type _bfield;
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
    // using surfaces_finder_t = typename detector_type::surfaces_finder_type;
    using bfield_t = typename detector_type::bfield_type::view_t;

    detector_data(detector_type &det)
        : _volumes_data(vecmem::get_data(det.volumes())),
          _surfaces_data(vecmem::get_data(det.surfaces())),
          _masks_data(get_data(det.mask_store())),
          _materials_data(get_data(det.material_store())),
          _transforms_data(get_data(det.transform_store())),
          _sf_finder_data(get_data(det.sf_finder_store())),
          _volume_finder_data(
              get_data(det.volume_search_grid(), *det.resource())),
          _bfield_data(det.get_bfield()) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    typename detector_type::mask_container::view_type _masks_data;
    typename detector_type::material_container::view_type _materials_data;
    typename detector_type::transform_container::view_type _transforms_data;
    typename detector_type::sf_finder_container::view_type _sf_finder_data;
    grid2_data<volume_finder_t> _volume_finder_data;
    bfield_t _bfield_data;
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
    using bfield_t = typename detector_type::bfield_type::view_t;

    detector_view(detector_data<detector_type> &det_data)
        : _volumes_data(det_data._volumes_data),
          _surfaces_data(det_data._surfaces_data),
          _masks_data(det_data._masks_data),
          _materials_data(det_data._materials_data),
          _transforms_data(det_data._transforms_data),
          _sf_finder_data(det_data._sf_finder_data),
          _volume_finder_view(det_data._volume_finder_data),
          _bfield_view(det_data._bfield_data) {}

    // members
    vecmem::data::vector_view<volume_t> _volumes_data;
    vecmem::data::vector_view<surface_t> _surfaces_data;
    typename detector_type::mask_container::view_type _masks_data;
    typename detector_type::material_container::view_type _materials_data;
    typename detector_type::transform_container::view_type _transforms_data;
    typename detector_type::sf_finder_container::view_type _sf_finder_data;
    grid2_view<volume_finder_t> _volume_finder_view;
    bfield_t _bfield_view;
};

/// stand alone function for that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename detector_registry, template <typename> class bfield_t,
          typename container_t, typename source_link>
inline detector_data<
    detector<detector_registry, bfield_t, container_t, source_link> >
get_data(detector<detector_registry, bfield_t, container_t, source_link> &det) {
    return det;
}

}  // namespace detray
