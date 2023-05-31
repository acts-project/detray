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
#include "detray/core/detail/detector_kernel.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
#include <algorithm>
#include <map>
#include <sstream>
#include <string>

namespace detray {

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
    using surface_container =
        typename metadata::template surface_finder_store<tuple_type,
                                                         container_t>;
    using sf_finders = typename surface_container::value_types;
    using sf_finder_link = typename surface_container::single_link;

    // TODO: Move to the following to volume builder

    /// The surface takes a mask (defines the local coordinates and the surface
    /// extent), its material, a link to an element in the transform container
    /// to define its placement and a source link to the object it represents.
    using surface_type = typename metadata::surface_type;

    using surface_container_t = vector_type<surface_type>;
    /// Volume type
    using geo_obj_ids = typename metadata::geo_objects;
    using volume_type =
        detector_volume<geo_obj_ids, sf_finder_link, scalar_type>;
    using volume_container = vector_type<volume_type>;

    /// Volume finder definition: Make volume index available from track
    /// position
    using volume_finder =
        typename metadata::template volume_finder<container_t>;

    using detector_view_type =
        detector_view<metadata, covfie::field, host_container_types>;

    using free_vector_type =
        typename free_track_parameters<transform3>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<transform3>::vector_type;

    /// Detector view types
    using view_type = dmulti_view<dvector_view<volume_type>,
                                  typename transform_container::view_type,
                                  typename mask_container::view_type,
                                  typename material_container::view_type,
                                  typename surface_container::view_type,
                                  typename volume_finder::view_type>;

    using const_view_type =
        dmulti_view<dvector_view<const volume_type>,
                    typename transform_container::const_view_type,
                    typename mask_container::const_view_type,
                    typename material_container::const_view_type,
                    typename surface_container::const_view_type,
                    typename volume_finder::const_view_type>;

    /// Detector buffer types
    using buffer_type = dmulti_buffer<dvector_buffer<volume_type>,
                                      typename transform_container::buffer_type,
                                      typename mask_container::buffer_type,
                                      typename material_container::buffer_type,
                                      typename surface_container::buffer_type,
                                      typename volume_finder::buffer_type>;

    detector() = delete;

    /// Allowed costructor
    /// @param resource memory resource for the allocation of members
    DETRAY_HOST
    detector(vecmem::memory_resource &resource, bfield_type &&field)
        : _volumes(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _surfaces(resource),
          _volume_finder(resource),
          _resource(&resource),
          _bfield(field) {}

    /// Constructor with simplified constant-zero B-field
    /// @param resource memory resource for the allocation of members
    DETRAY_HOST
    explicit detector(vecmem::memory_resource &resource)
        : _volumes(&resource),
          _transforms(resource),
          _masks(resource),
          _materials(resource),
          _surfaces(resource),
          _volume_finder(resource),
          _resource(&resource),
          _bfield(typename bfield_type::backend_t::configuration_t{0.f, 0.f,
                                                                   0.f}) {}

    /// Constructor with detector_data
    template <typename detector_data_type,
              std::enable_if_t<!std::is_base_of_v<vecmem::memory_resource,
                                                  detector_data_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE explicit detector(detector_data_type &det_data)
        : _volumes(detray::detail::get<0>(det_data._detector_data.m_view)),
          _transforms(detray::detail::get<1>(det_data._detector_data.m_view)),
          _masks(detray::detail::get<2>(det_data._detector_data.m_view)),
          _materials(detray::detail::get<3>(det_data._detector_data.m_view)),
          _surfaces(detray::detail::get<4>(det_data._detector_data.m_view)),
          _volume_finder(
              detray::detail::get<5>(det_data._detector_data.m_view)),
          _bfield(det_data._bfield_view) {}

    /// Add a new volume and retrieve a reference to it.
    ///
    /// @param id the shape id for the volume
    /// @param sf_finder_link of the volume, where to entry the surface finder
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

    /// Add a new volume and retrieve a reference to it.
    ///
    /// @param id the shape id for the volume
    /// @param bounds of the volume, they are expected to be already attaching
    /// @param sf_finder_link of the volume, where to entry the surface finder
    ///
    /// @return non-const reference to the new volume
    DETRAY_HOST
    volume_type &new_volume(
        const volume_id id, const array_type<scalar, 6> &bounds,
        typename volume_type::link_type::index_type srf_finder_link = {}) {
        volume_type &cvolume = new_volume(id, srf_finder_link);
        cvolume.set_bounds(bounds);

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

    /// @return the volume by global cartesian @param position - const access
    DETRAY_HOST_DEVICE
    inline auto volume_by_pos(const point3 &p) const -> const volume_type & {
        // The 3D cylindrical volume search grid is concentric
        const transform3 identity{};
        const auto loc_pos =
            _volume_finder.global_to_local(identity, p, identity.translation());

        // Only one entry per bin
        dindex volume_index{*_volume_finder.search(loc_pos)};
        return _volumes[volume_index];
    }

    /// @returns access to the surface finder container
    DETRAY_HOST_DEVICE
    inline auto surface_store() const -> const surface_container & {
        return _surfaces;
    }

    /// @returns access to the surface finder container
    DETRAY_HOST_DEVICE
    inline auto surface_store() -> surface_container & { return _surfaces; }

    /// @returns all surfaces - const
    DETRAY_HOST_DEVICE
    inline const auto &surfaces() const {
        return _surfaces.template get<sf_finders::id::e_brute_force>().all();
    }

    /// @returns all surfaces - non-const
    DETRAY_HOST_DEVICE
    inline auto &surfaces() {
        return _surfaces.template get<sf_finders::id::e_brute_force>().all();
    }

    /// @returns a surface by index - const
    DETRAY_HOST_DEVICE
    constexpr auto surfaces(geometry::barcode bcd) const
        -> const surface_type & {
        return _surfaces.template get<sf_finders::id::e_brute_force>()
            .all()[bcd.index()];
    }

    /// @returns surfaces of a given type (@tparam sf_id) by volume - const
    template <geo_obj_ids sf_id = static_cast<geo_obj_ids>(0)>
    DETRAY_HOST_DEVICE constexpr auto surfaces(const volume_type &v) const {
        dindex coll_idx{v.template link<sf_id>().index()};
        const auto sf_coll =
            _surfaces.template get<sf_finders::id::e_brute_force>()[coll_idx];
        return sf_coll.all();
    }

    /// Append new surfaces to the detector
    DETRAY_HOST
    inline void append_surfaces(surface_container_t &&new_surfaces) {
        /*_surfaces.insert(_surfaces.end(), new_surfaces.begin(),
                         new_surfaces.end());*/
        _surfaces.template push_back<sf_finders::id::e_brute_force>(
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
        surface_container_t &surfaces_per_vol, mask_container &masks_per_vol,
        transform_container &trfs_per_vol) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = _transforms.size(ctx);
        _transforms.append(std::move(trfs_per_vol), ctx);

        // Update mask, material and transform index of surfaces and set a
        // unique barcode (index of surface in container)
        dindex sf_offset{static_cast<dindex>(
            _surfaces.template get<sf_finders::id::e_brute_force>()
                .all()
                .size())};
        for (auto &sf : surfaces_per_vol) {
            _masks.template visit<detail::mask_index_update>(sf.mask(), sf);
            sf.update_transform(trf_offset);
            sf.set_index(sf_offset++);
        }

        // Append surfaces to base surface collection
        _surfaces.template push_back<sf_finders::id::e_brute_force>(
            surfaces_per_vol);

        // Update the surface link in a volume
        vol.template set_link<surface_id>(
            sf_finders::id::e_brute_force,
            _surfaces.template size<sf_finders::id::e_brute_force>() - 1);

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
        surface_container_t &surfaces_per_vol, mask_container &masks_per_vol,
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

    /// @returns the maximum number of surface candidates that any volume may
    /// return.
    DETRAY_HOST
    inline auto n_max_candidates() const -> std::size_t {
        std::vector<std::size_t> n_candidates;
        n_candidates.reserve(_volumes.size());
        for (const auto &vol : _volumes) {
            n_candidates.push_back(vol.n_max_candidates(surface_store()));
        }
        return *std::max_element(n_candidates.begin(), n_candidates.end());
    }

    DETRAY_HOST_DEVICE
    inline const bfield_type &get_bfield() const { return _bfield; }

    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const geometry::barcode bc, const point3 &pos,
                                  const vector3 &dir) const {
        const auto &sf =
            *(surfaces().begin() + static_cast<std::ptrdiff_t>(bc.index()));
        const auto ret =
            _masks.template visit<detail::global_to_local<transform3>>(
                sf.mask(), _transforms[sf.transform()], pos, dir);
        return ret;
    }

    DETRAY_HOST_DEVICE
    inline point3 local_to_global(const geometry::barcode bc,
                                  const point3 &local) const {
        const auto &sf =
            *(surfaces().begin() + static_cast<std::ptrdiff_t>(bc.index()));
        const auto ret =
            _masks.template visit<detail::local_to_global<transform3>>(
                sf.mask(), _transforms[sf.transform()], local);
        return ret;
    }

    DETRAY_HOST_DEVICE
    inline bound_vector_type free_to_bound_vector(
        const geometry::barcode bc, const free_vector_type &free_vec) const {
        const auto &sf =
            *(surfaces().begin() + static_cast<std::ptrdiff_t>(bc.index()));
        const auto ret =
            _masks.template visit<detail::free_to_bound_vector<transform3>>(
                sf.mask(), _transforms[sf.transform()], free_vec);
        return ret;
    }

    DETRAY_HOST_DEVICE
    inline free_vector_type bound_to_free_vector(
        const geometry::barcode bc, const bound_vector_type &bound_vec) const {
        const auto &sf =
            *(surfaces().begin() + static_cast<std::ptrdiff_t>(bc.index()));
        const auto ret =
            _masks.template visit<detail::bound_to_free_vector<transform3>>(
                sf.mask(), _transforms[sf.transform()], bound_vec);
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

            ss << "                 " << _surfaces.n_collections()
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
    volume_container _volumes;

    /// Keeps all of the transform data in contiguous memory
    transform_container _transforms;

    /// Masks of all surfaces in the geometry in contiguous memory
    mask_container _masks;

    /// Materials of all surfaces in the geometry in contiguous memory
    material_container _materials;

    /// All surface finder data structures that are used in the detector volumes
    surface_container _surfaces;

    /// Search structure for volumes
    volume_finder _volume_finder;

    /// The memory resource represents how and where (host, device, managed)
    /// the memory for the detector containers is allocated
    vecmem::memory_resource *_resource = nullptr;

    /// Storage for magnetic field data
    bfield_type _bfield;
};

/// @brief A static inplementation of detector data for device
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
struct detector_buffer {

    using detector_type = detector<metadata, bfield_t, container_t>;

    /// Automatic buffer creation with the given parameters
    detector_buffer(detector_type &det, vecmem::memory_resource &mr,
                    vecmem::copy &cpy,
                    detray::copy cpy_type = detray::copy::sync,
                    vecmem::data::buffer_type buff_type =
                        vecmem::data::buffer_type::fixed_size)
        : _detector_buffer(detray::get_buffer(detray::get_data(det.volumes()),
                                              mr, cpy, cpy_type, buff_type),
                           detray::get_buffer(det.transform_store(), mr, cpy,
                                              cpy_type, buff_type),
                           detray::get_buffer(det.mask_store(), mr, cpy,
                                              cpy_type, buff_type),
                           detray::get_buffer(det.material_store(), mr, cpy,
                                              cpy_type, buff_type),
                           detray::get_buffer(det.surface_store(), mr, cpy,
                                              cpy_type, buff_type),
                           detray::get_buffer(det.volume_search_grid(), mr, cpy,
                                              cpy_type, buff_type)),
          _bfield_view(det.get_bfield()) {}

    /// Buffers were created manually
    detector_buffer(
        detector<metadata, bfield_t, container_t> &det,
        detail::get_buffer_t<typename detector_type::volume_container>
            &&vol_buffer,
        typename detector_type::transform_container::buffer_type &&trf_buffer,
        typename detector_type::mask_container::buffer_type &&msk_buffer,
        typename detector_type::material_container::buffer_type &&mat_buffer,
        typename detector_type::surface_container::buffer_type &&sf_buffer,
        typename detector_type::volume_finder::buffer_type &&vgrd_buffer)
        : _detector_buffer(std::move(vol_buffer), std::move(trf_buffer),
                           std::move(msk_buffer), std::move(mat_buffer),
                           std::move(sf_buffer), std::move(vgrd_buffer)),
          _bfield_view(det.get_bfield()) {}

    /// Buffers for the vecemem types
    typename detector_type::buffer_type _detector_buffer;
    /// Covfie field
    typename detector_type::bfield_type::view_t _bfield_view;
};

// Deduction guide to construct a detector buffer type more conveniently
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
detector_buffer(
    detector<metadata, bfield_t, container_t> &,
    detail::get_buffer_t<typename detector<metadata, bfield_t,
                                           container_t>::volume_container> &&,
    typename detector<metadata, bfield_t,
                      container_t>::transform_container::buffer_type &&,
    typename detector<metadata, bfield_t,
                      container_t>::mask_container::buffer_type &&,
    typename detector<metadata, bfield_t,
                      container_t>::material_container::buffer_type &&,
    typename detector<metadata, bfield_t,
                      container_t>::surface_container::buffer_type &&,
    typename detector<metadata, bfield_t,
                      container_t>::volume_finder::buffer_type &&)
    -> detector_buffer<metadata, bfield_t, container_t>;

/// @brief A static inplementation of detector data for device
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
struct detector_view {

    using detector_type = detector<metadata, bfield_t, container_t>;

    detector_view(detector_type &det)
        : _detector_data(detray::get_data(det.volumes()),
                         detray::get_data(det.transform_store()),
                         detray::get_data(det.mask_store()),
                         detray::get_data(det.material_store()),
                         detray::get_data(det.surface_store()),
                         detray::get_data(det.volume_search_grid())),
          _bfield_view(det.get_bfield()) {}

    detector_view(detector_buffer<metadata, bfield_t, container_t> &det_buff)
        : _detector_data(detray::get_data(det_buff._detector_buffer)),
          _bfield_view(det_buff._bfield_view) {}

    /// Views for the vecmem types
    typename detector_type::view_type _detector_data;
    /// Covfie field view
    typename detector_type::bfield_type::view_t _bfield_view;
};

/// Stand-alone function that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
inline detector_view<metadata, bfield_t, container_t> get_data(
    detector<metadata, bfield_t, container_t> &det) {
    return {det};
}

/// Stand-alone function that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
inline detector_buffer<metadata, bfield_t, container_t> get_buffer(
    detector<metadata, bfield_t, container_t> &det, vecmem::memory_resource &mr,
    vecmem::copy &cpy, detray::copy cpy_type = detray::copy::sync,
    vecmem::data::buffer_type buff_type =
        vecmem::data::buffer_type::fixed_size) {
    return {det, mr, cpy, cpy_type, buff_type};
}

/// Stand-alone function that @returns the detector data for transfer to
/// device.
///
/// @param detector the detector to be tranferred
template <typename metadata, template <typename> class bfield_t,
          typename container_t>
inline detector_view<metadata, bfield_t, container_t> get_data(
    detector_buffer<metadata, bfield_t, container_t> &det_buff) {
    return {det_buff};
}

}  // namespace detray
