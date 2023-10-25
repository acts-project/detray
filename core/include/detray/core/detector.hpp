/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detail/volume_descriptor.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <algorithm>
#include <map>
#include <sstream>
#include <string>

namespace detray {

/// @brief Forward declaration of a detector view type
template <typename metadata>
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
template <typename metadata_t = default_metadata,
          typename container_t = host_container_types>
class detector {

    // Allow the building of the detector containers
    friend class volume_builder_interface<detector<metadata_t, container_t>>;

    public:
    /// Algebra types
    /// @TODO: scalar as a template parameter
    using scalar_type = scalar;

    using point3 = __plugin::point3<scalar_type>;
    using vector3 = __plugin::vector3<scalar_type>;
    using point2 = __plugin::point2<scalar_type>;

    using metadata = metadata_t;

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

    /// The surface takes a mask (defines the local coordinates and the surface
    /// extent), its material, a link to an element in the transform container
    /// to define its placement and a source link to the object it represents.
    using surface_type = typename metadata::surface_type;
    using surface_container = vector_type<surface_type>;
    using surface_lookup_container = surface_container;

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
    using accelerator_container =
        typename metadata::template accelerator_store<tuple_type, container_t>;
    using accel = typename accelerator_container::value_types;
    using accel_link = typename accelerator_container::single_link;

    /// Volume type
    using geo_obj_ids = typename metadata::geo_objects;
    using volume_type = volume_descriptor<geo_obj_ids, accel_link>;
    using volume_container = vector_type<volume_type>;

    /// Volume finder definition: Make volume index available from track
    /// position
    using volume_finder =
        typename metadata::template volume_finder<container_t>;

    /// Detector view types
    using view_type =
        dmulti_view<dvector_view<volume_type>, dvector_view<surface_type>,
                    typename transform_container::view_type,
                    typename mask_container::view_type,
                    typename material_container::view_type,
                    typename accelerator_container::view_type,
                    typename volume_finder::view_type>;

    static_assert(detail::is_device_view_v<view_type>,
                  "Detector view type ill-formed");

    using const_view_type =
        dmulti_view<dvector_view<const volume_type>,
                    dvector_view<const surface_type>,
                    typename transform_container::const_view_type,
                    typename mask_container::const_view_type,
                    typename material_container::const_view_type,
                    typename accelerator_container::const_view_type,
                    typename volume_finder::const_view_type>;

    static_assert(detail::is_device_view_v<const_view_type>,
                  "Detector const view type ill-formed");

    /// Detector buffer types
    using buffer_type =
        dmulti_buffer<dvector_buffer<volume_type>, dvector_buffer<surface_type>,
                      typename transform_container::buffer_type,
                      typename mask_container::buffer_type,
                      typename material_container::buffer_type,
                      typename accelerator_container::buffer_type,
                      typename volume_finder::buffer_type>;

    static_assert(detail::is_buffer_v<buffer_type>,
                  "Detector buffer type ill-formed");

    detector() = delete;
    // The detector holds a lot of data and should never be copied
    detector(const detector &) = delete;
    detector &operator=(const detector &) = delete;
    detector(detector &&) = default;

    /// Allowed costructor
    /// @{
    /// Default construction
    /// @param resource memory resource for the allocation of members
    DETRAY_HOST
    explicit detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _surfaces(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _accelerators(resource),
          _volume_finder(resource),
          _resource(&resource) {}

    /// Constructor with detector_data
    template <typename detector_view_t,
              typename std::enable_if_t<
                  detail::is_device_view_v<detector_view_t>, bool> = true>
    DETRAY_HOST_DEVICE explicit detector(detector_view_t &det_data)
        : _volumes(detray::detail::get<0>(det_data.m_view)),
          _surfaces(detray::detail::get<1>(det_data.m_view)),
          _transforms(detray::detail::get<2>(det_data.m_view)),
          _masks(detray::detail::get<3>(det_data.m_view)),
          _materials(detray::detail::get<4>(det_data.m_view)),
          _accelerators(detray::detail::get<5>(det_data.m_view)),
          _volume_finder(detray::detail::get<6>(det_data.m_view)) {}
    /// @}

    /// Add a new volume and retrieve a reference to it.
    ///
    /// @param id the shape id for the volume
    /// @param accel_link of the volume, where to entry the surface finder
    ///
    /// @return non-const reference to the new volume
    DETRAY_HOST
    volume_type &new_volume(
        const volume_id id,
        typename volume_type::link_type::index_type srf_finder_link = {}) {
        volume_type &cvolume = _volumes.emplace_back(id);
        cvolume.set_index(static_cast<dindex>(_volumes.size()) - 1u);
        cvolume
            .template set_link<static_cast<typename volume_type::object_id>(0)>(
                srf_finder_link);

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

    /// @return the volume for a @param volume_descr - const access
    DETRAY_HOST_DEVICE
    inline auto volume(const volume_type &volume_descr) const
        -> const detector_volume<detector> {
        return detector_volume{*this, volume_descr};
    }

    /// @return the volume for a @param volume_descr - non-const access
    DETRAY_HOST_DEVICE
    inline auto volume(volume_type &volume_descr) -> detector_volume<detector> {
        return detector_volume{*this, volume_descr};
    }

    /// @return the volume by @param volume_index - const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_index(dindex volume_index) const
        -> const detector_volume<detector> {
        return detector_volume{*this, _volumes[volume_index]};
    }

    /// @return the volume by @param volume_index - non-const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_index(dindex volume_index)
        -> detector_volume<detector> {
        return detector_volume{*this, _volumes[volume_index]};
    }

    /// @return the volume by global cartesian @param position - const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_pos(const point3 &p) const
        -> const detector_volume<detector> {
        // The 3D cylindrical volume search grid is concentric
        const transform3 identity{};
        const auto loc_pos =
            _volume_finder.global_to_local(identity, p, identity.translation());

        // Only one entry per bin
        dindex volume_index{*_volume_finder.search(loc_pos)};
        return detector_volume{*this, _volumes[volume_index]};
    }

    /// @returns access to the surface finder container
    DETRAY_HOST_DEVICE
    inline auto accelerator_store() const -> const accelerator_container & {
        return _accelerators;
    }

    /// @returns access to the surface finder container
    DETRAY_HOST_DEVICE
    inline auto accelerator_store() -> accelerator_container & {
        return _accelerators;
    }

    /// @returns all portals - const
    /// @note Depending on the detector type, this can also contain other
    /// surfaces
    DETRAY_HOST_DEVICE
    inline const auto &portals() const {
        // In case of portals, we know where they live
        return _accelerators.template get<accel::id::e_brute_force>().all();
    }

    /// @returns all portals - non-const
    /// @note Depending on the detector type, this can also container other
    /// surfaces
    DETRAY_HOST_DEVICE
    inline auto &portals() {
        return _accelerators.template get<accel::id::e_brute_force>().all();
    }

    /// @returns the portals of a given volume @param v - const
    /// @note Depending on the detector type, this can also container other
    /// surfaces
    DETRAY_HOST_DEVICE constexpr auto portals(const volume_type &v) const {

        // Index of the portal collection in the surface store
        const dindex pt_coll_idx{
            v.template link<geo_obj_ids::e_portal>().index()};

        const auto &pt_coll =
            _accelerators.template get<accel::id::e_brute_force>()[pt_coll_idx];

        return pt_coll.all();
    }

    /// @return the sub-volumes of the detector - const access
    DETRAY_HOST_DEVICE
    inline auto surface_lookup() const -> const vector_type<surface_type> & {
        return _surfaces;
    }

    /// @return the sub-volumes of the detector - non-const access
    DETRAY_HOST_DEVICE
    inline auto surface_lookup() -> vector_type<surface_type> & {
        return _surfaces;
    }

    /// @returns a surface using its barcode - const
    DETRAY_HOST_DEVICE
    constexpr auto surface(geometry::barcode bcd) const
        -> const surface_type & {
        return _surfaces[bcd.index()];
    }

    /// @returns the overall number of surfaces in the detector
    DETRAY_HOST_DEVICE
    constexpr auto n_surfaces() const -> dindex {
        return static_cast<dindex>(_surfaces.size());
    }

    /// Add a new surface to the lookup according to its index.
    DETRAY_HOST
    constexpr auto add_surface_to_lookup(const surface_type sf) -> void {
        if (_surfaces.size() < sf.index() + 1) {
            _surfaces.resize(sf.index() + 1);
        }
        _surfaces.at(sf.index()) = sf;
    }

    /// Append new portals(surfaces) to the detector
    DETRAY_HOST
    inline void append_portals(surface_container &&new_surfaces) {
        _accelerators.template push_back<accel::id::e_brute_force>(
            std::move(new_surfaces));
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
    template <geo_obj_ids surface_id = static_cast<geo_obj_ids>(0)>
    DETRAY_HOST auto add_objects_per_volume(
        const geometry_context ctx, volume_type &vol,
        surface_container &surfaces_per_vol, mask_container &masks_per_vol,
        transform_container &trfs_per_vol) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = _transforms.size(ctx);
        _transforms.append(std::move(trfs_per_vol), ctx);

        // Update mask, material and transform index of surfaces and set a
        // unique barcode (index of surface in container)
        for (auto &sf : surfaces_per_vol) {
            _masks.template visit<detail::mask_index_update>(sf.mask(), sf);
            sf.update_transform(trf_offset);
            // Don't overwrite the surface index if it has been set already
            // (e.g. with a special sorting)
            if (sf.barcode().is_invalid()) {
                sf.set_index(n_surfaces());
                assert(!sf.barcode().is_invalid());
            }
            add_surface_to_lookup(sf);
        }

        // Append surfaces to base surface collection
        _accelerators.template push_back<accel::id::e_default>(
            surfaces_per_vol);

        // Update the surface link in a volume
        vol.template set_link<surface_id>(
            accel::id::e_default,
            _accelerators.template size<accel::id::e_default>() - 1);

        // Append mask and material container
        _masks.append(std::move(masks_per_vol));
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
            _materials.template visit<detail::material_index_update>(
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
    inline auto set_volume_finder(volume_finder &&v_grid) -> void {
        _volume_finder = std::move(v_grid);
    }

    /// Add the volume grid - copy semantics
    ///
    /// @param v_grid the volume grid to be added
    DETRAY_HOST
    inline auto set_volume_finder(const volume_finder &v_grid) -> void {
        _volume_finder = v_grid;
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

    /// @returns the maximum number of surface candidates that any volume may
    /// return.
    DETRAY_HOST
    inline auto n_max_candidates() const -> std::size_t {
        std::vector<std::size_t> n_candidates;
        n_candidates.reserve(_volumes.size());
        for (const auto &vol : _volumes) {
            n_candidates.push_back(
                detector_volume{*this, vol}.n_max_candidates());
        }
        return *std::max_element(n_candidates.begin(), n_candidates.end());
    }

    /// @returns view of a detector
    DETRAY_HOST auto get_data() -> view_type {
        return view_type{
            detray::get_data(_volumes),      detray::get_data(_surfaces),
            detray::get_data(_transforms),   detray::get_data(_masks),
            detray::get_data(_materials),    detray::get_data(_accelerators),
            detray::get_data(_volume_finder)};
    }

    /// @returns const view of a detector
    DETRAY_HOST auto get_data() const -> const_view_type {
        return const_view_type{
            detray::get_data(_volumes),      detray::get_data(_surfaces),
            detray::get_data(_transforms),   detray::get_data(_masks),
            detray::get_data(_materials),    detray::get_data(_accelerators),
            detray::get_data(_volume_finder)};
    }

    /// @param names maps a volume to its string representation.
    /// @returns a string representation of the detector.
    DETRAY_HOST
    auto to_string(const name_map &names) const -> std::string {
        std::stringstream ss;

        ss << "[>] Detector '" << names.at(0) << "' has " << _volumes.size()
           << " volumes." << std::endl;

        for (const auto [i, v] : detray::views::enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << ": " << std::endl;
            ss << " - name: '" << names.at(v.index() + 1u) << "'" << std::endl;

            ss << "     contains    "
               << v.template n_objects<geo_obj_ids::e_sensitive>()
               << " sensitive surfaces " << std::endl;

            ss << "                 "
               << v.template n_objects<geo_obj_ids::e_portal>() << " portals "
               << std::endl;

            ss << "                 " << _accelerators.n_collections()
               << " surface finders " << std::endl;

            if (v.accel_index() != dindex_invalid) {
                ss << "  sf finder id " << v.accel_type() << "  sf finders idx "
                   << v.accel_index() << std::endl;
            }
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
    volume_container _volumes;

    /// Lookup for surfaces from barcodes
    surface_lookup_container _surfaces;

    /// Keeps all of the transform data in contiguous memory
    transform_container _transforms;

    /// Masks of all surfaces in the geometry in contiguous memory
    mask_container _masks;

    /// Materials of all surfaces in the geometry in contiguous memory
    material_container _materials;

    /// All surface finder data structures that are used in the detector volumes
    accelerator_container _accelerators;

    /// Search structure for volumes
    volume_finder _volume_finder;

    /// The memory resource represents how and where (host, device, managed)
    /// the memory for the detector containers is allocated
    vecmem::memory_resource *_resource = nullptr;
};

}  // namespace detray
