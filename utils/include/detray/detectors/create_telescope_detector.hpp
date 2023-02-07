/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/line_stepper.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
#include <limits>

namespace detray {

namespace {

// @Note: These type definitions should be removed at some point
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

template <typename mask_shape_t>
using telescope_types =
    typename detector_registry::template telescope_detector<mask_shape_t>;

/// Where and how to place the telescope modules.
struct module_placement {
    /// Module position
    point3 _pos;
    /// Module normal
    vector3 _dir;

    bool operator==(const module_placement &other) const {
        bool is_same = true;
        is_same &= (_pos == other._pos);
        is_same &= (_dir == other._dir);
        return is_same;
    }
};

/// Helper method for positioning the plane surfaces.
///
/// @param traj pilot trajectory along which the modules should be placed.
/// @param steps lengths along the trajectory where surfaces should be placed.
///
/// @return a vector of the @c module_placements along the trajectory.
template <typename trajectory_t>
inline std::vector<module_placement> module_positions(
    const trajectory_t &traj, const std::vector<scalar> &steps) {

    // create and fill the module placements
    std::vector<module_placement> placements;
    placements.reserve(steps.size());

    for (const auto s : steps) {
        placements.push_back({traj.pos(s), traj.dir(s)});
    }

    return placements;
}

/// Helper function that creates the telescope surfaces.
///
/// @param ctx geometric context
/// @param traj pilot trajectory along which the modules should be placed
/// @param volume volume the planes should be added to
/// @param surfaces container to add new surface to
/// @param masks container to add new cylinder mask to
/// @param transforms container to add new transform to
/// @param cfg config struct for module creation
template <typename mask_t, typename context_t, typename trajectory_t,
          typename volume_type, typename surface_container_t,
          typename mask_container_t, typename material_container_t,
          typename transform_container_t, typename config_t>
inline void create_telescope(context_t &ctx, const trajectory_t &traj,
                             volume_type &volume, surface_container_t &surfaces,
                             mask_container_t &masks,
                             material_container_t &materials,
                             transform_container_t &transforms,
                             const config_t &cfg) {
    using surface_type = typename surface_container_t::value_type;
    using volume_link_t = typename surface_type::volume_link_type;
    using mask_link_type = typename surface_type::mask_link;
    using material_link_type = typename surface_type::material_link;

    auto volume_idx = volume.index();
    constexpr auto slab_id = material_link_type::id_type::e_slab;
    constexpr typename mask_container_t::ids mask_id{0};

    // Create the module centers
    const std::vector<module_placement> m_placements =
        module_positions(traj, cfg.pos);

    // Create geometry data
    for (const auto &m_placement : m_placements) {

        volume_link_t mask_volume_link{volume_idx};

        // Surfaces with the linking into the local containers
        mask_link_type mask_link{mask_id, masks.template size<mask_id>()};
        material_link_type material_link{slab_id,
                                         materials.template size<slab_id>()};
        const auto trf_index = transforms.size(ctx);
        surfaces.emplace_back(trf_index, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_sensitive);

        // The first and last surface acts as portal that leaves the telescope
        if (m_placement == m_placements.front()) {
            mask_volume_link = dindex_invalid;
            surfaces.front().set_id(surface_id::e_portal);
        }

        if (m_placement == m_placements.back()) {
            mask_volume_link = dindex_invalid;
            surfaces.back().set_id(surface_id::e_portal);
        }

        // The rectangle bounds for this module
        masks.template emplace_back<mask_id>(empty_context{}, cfg.mask_values,
                                             mask_volume_link);

        // This will not catch every possible line type...
        if constexpr (std::is_same_v<typename mask_t::shape, line<>> or
                      std::is_same_v<typename mask_t::shape, line<true>>) {
            materials.template emplace_back<material_container_t::ids::e_rod>(
                empty_context{}, cfg.m_mat, cfg.m_thickness);
        } else {
            materials.template emplace_back<material_container_t::ids::e_slab>(
                empty_context{}, cfg.m_mat, cfg.m_thickness);
        }

        // Build the transform
        // Local z axis is the global normal vector
        vector3 m_local_z = algebra::vector::normalize(m_placement._dir);

        // Project onto the weakest direction component of the normal vector
        vector3 e_i{0.f, 0.f, 0.f};
        scalar min = std::numeric_limits<scalar>::infinity();
        unsigned int i{std::numeric_limits<uint>::max()};
        for (unsigned int k = 0u; k < 3u; ++k) {
            if (m_local_z[k] < min) {
                min = m_local_z[k];
                i = k;
            }
        }
        e_i[i] = 1.f;
        vector3 proj = algebra::vector::dot(m_local_z, e_i) * m_local_z;
        // Local x axis is the normal to local y,z
        vector3 m_local_x = algebra::vector::normalize(e_i - proj);

        // Create the global-to-local transform of the module
        transforms.emplace_back(ctx, m_placement._pos, m_local_z, m_local_x);
    }
}

}  // namespace

/// Builds a detray geometry that contains only one volume with one type of
/// surfaces, where the last surface is the portal that leaves the telescope.
/// The detector is auto-constructed by following a trajectory state through
/// space. The trajectory and the given distances determine the positions of
/// the plane surfaces. The dis
///
/// @tparam mask_t the type of mask for the telescope surfaces
/// @tparam trajectory_t the type of the pilot trajectory
/// @tparam container_t the containers used in the detector data structure
/// @tparam Args value type for the mask boundaries
///
/// @param resource the memory resource for the detector containers
/// @param bfield the magnetic field description for the detector
/// @param pos the module positions. These only correspond to the actual module
///            positions for a straight line pilot trajectory
/// @param traj the pilot trajectory along which the surfaces are positioned
/// @param mat the surface material
/// @param thickness the thisckness of the surface material (slab/rod)
/// @param msk the surface boundary mask
///
/// @returns a complete detector object
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types>
auto create_telescope_detector(
    vecmem::memory_resource &resource,
<<<<<<< HEAD
    covfie::field<detector_registry::telescope_detector::bfield_backend_t>
        &&bfield,
    std::vector<scalar> pos,
    trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    scalar half_x = 20.f * unit<scalar>::mm,
    scalar half_y = 20.f * unit<scalar>::mm,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um) {
=======
    covfie::field<typename detector_registry::template telescope_detector<
        typename mask_t::shape>::bfield_backend_t> &&bfield,
    const mask_t &msk, std::vector<scalar> pos,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80 * unit<scalar>::um,
    trajectory_t traj = {{0, 0, 0}, 0, {0, 0, 1}, -1}) {
>>>>>>> a1f3d421 (adapt telescope detector to hold single but different mask types)

    // detector type
    using detector_t = detector<telescope_types<typename mask_t::shape>,
                                covfie::field, container_t>;

    // module parameters
    struct surface_config {
        const typename mask_t::mask_values &mask_values;
        const std::vector<scalar> &pos;
        material<scalar> m_mat;
        scalar m_thickness;
    };

    // create empty detector
    detector_t det(resource, std::move(bfield));

    typename detector_t::geometry_context ctx{};
    const surface_config sf_config{msk.values(), pos, mat, thickness};

    // volume boundaries are not needed. Same goes for portals
    det.new_volume(
        volume_id::e_cylinder,
        {0.f, 0.f, 0.f, 0.f, -constant<scalar>::pi, constant<scalar>::pi});
    typename detector_t::volume_type &vol = det.volume_by_index(0u);

    // Add module surfaces to volume
    typename detector_t::surface_container surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    create_telescope<mask_t>(ctx, traj, vol, surfaces, masks, materials,
                             transforms, sf_config);

    det.add_objects_per_volume(ctx, vol, surfaces, masks, transforms,
                               materials);

    return det;
}

/// Build the telescope geometry from a fixed length and number of surfaces.
///
/// @param n_surfaces the number of surfaces that are placed in the geometry
/// @param tel_length the total length of the steps by the stepper
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types, typename... Args>
auto create_telescope_detector(
    vecmem::memory_resource &resource,
    covfie::field<
        typename telescope_types<typename mask_t::shape>::bfield_backend_t>
        &&bfield,
<<<<<<< HEAD
    std::size_t n_surfaces = 10u, scalar tel_length = 500.f * unit<scalar>::mm,
    trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    scalar half_x = 20.f * unit<scalar>::mm,
    scalar half_y = 20.f * unit<scalar>::mm) {
=======
    const mask_t &msk, dindex n_surfaces = 10,
    scalar tel_length = 500. * unit<scalar>::mm,
    trajectory_t traj = {{0, 0, 0}, 0, {0, 0, 1}, -1}) {
>>>>>>> a1f3d421 (adapt telescope detector to hold single but different mask types)
    // Generate equidistant positions
    std::vector<scalar> positions = {};
    scalar pos{0.f};
    scalar dist{tel_length / static_cast<scalar>(n_surfaces - 1u)};
    for (std::size_t i = 0u; i < n_surfaces; ++i) {
        positions.push_back(pos);
        pos += dist;
    }

    // Set the default material
    const material<scalar> mat{silicon_tml<scalar>()};
    const scalar thickness{80.f * unit<scalar>::um};

    // Build the geometry
    return create_telescope_detector<mask_t>(resource, std::move(bfield), msk,
                                             positions, mat, thickness, traj);
}

/// Wrapper for create_telescope_geometry with constant zero bfield.
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types>
auto create_telescope_detector(
    vecmem::memory_resource &resource, const mask_t &msk,
    std::vector<scalar> pos, const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um,
    trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f}) {

    using covfie_bkdn_t =
        typename telescope_types<typename mask_t::shape>::bfield_backend_t;

    // Build the geometry
    return create_telescope_detector<mask_t, trajectory_t, container_t>(
        resource,
        covfie::field<covfie_bkdn_t>{
            typename covfie_bkdn_t::configuration_t{0.f, 0.f, 0.f}},
        msk, pos, mat, thickness, traj);
}

/// Wrapper for create_telescope_geometry with constant zero bfield.
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types>
auto create_telescope_detector(vecmem::memory_resource &resource,
                               const mask_t &msk, std::size_t n_surfaces = 10u,
                               scalar tel_length = 500.f * unit<scalar>::mm,
                               trajectory_t traj = {
                                   {0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f}) {

    using covfie_bkdn_t =
        typename telescope_types<typename mask_t::shape>::bfield_backend_t;

    // Build the geometry
    return create_telescope_detector<mask_t, trajectory_t, container_t>(
        resource,
        covfie::field<covfie_bkdn_t>{
            typename covfie_bkdn_t::configuration_t{0.f, 0.f, 0.f}},
        msk, n_surfaces, tel_length, traj);
}

}  // namespace detray
