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
#include "detray/masks/masks.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/tools/bounding_volume.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
#include <limits>

namespace detray {

namespace {

template <typename mask_shape_t>
using telescope_types =
    typename detector_registry::template telescope_detector<mask_shape_t>;

/// Where and how to place the telescope modules.
struct module_placement {

    // @Note: These type definitions should be removed at some point
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

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

/// Helper function to create minimum aabbs around all module surfaces in the
/// telescope and then construct world portals from the global bounding box.
/// This assumes that the module surfaces are already built inside the
/// argument containers.
///
/// @param ctx geometry context.
/// @param pt_envelope envelope around the module surfaces for the portals.
/// @param volume the single cuboid volume that the portals will be added to.
/// @param surfaces the module surface descriptors.
/// @param masks the module masks.
/// @param materials the materials (only needed to add the portal materials).
/// @param transforms the module surface transforms.
template <auto mask_id, typename context_t, typename volume_t,
          typename surface_container_t, typename mask_container_t,
          typename material_container_t, typename transform_container_t>
inline void create_cuboid_portals(context_t &ctx, const scalar pt_envelope,
                                  volume_t &volume,
                                  surface_container_t &surfaces,
                                  mask_container_t &masks,
                                  material_container_t &materials,
                                  transform_container_t &transforms) {
    // @Note: These type definitions should be removed at some point
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    using surface_t = typename surface_container_t::value_type;
    using nav_link_t = typename surface_t::navigation_link;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    using aabb_t = axis_aligned_bounding_volume<cuboid3D<>>;

    constexpr auto rectangle_id{mask_container_t::ids::e_portal_rectangle2};
    constexpr auto slab_id{material_container_t::ids::e_slab};

    // Envelope for the module surface minimum aabb
    constexpr scalar envelope{10.f * std::numeric_limits<scalar>::epsilon()};
    // Max distance in case of infinite bounds
    constexpr scalar max_shift{0.01f * std::numeric_limits<scalar>::max()};

    // The bounding boxes around the module surfaces
    std::vector<aabb_t> boxes;
    boxes.reserve(surfaces.size());

    const auto &mask_collection = masks.template get<mask_id>();

    for (const auto &sf : surfaces) {
        // Local minimum bounding box
        aabb_t box{mask_collection[sf.mask().index()], boxes.size(), envelope};
        // Minimum bounding box in global coordinates
        boxes.push_back(box.transform(transforms[sf.transform()]));
    }
    // Build an aabb in the global space around the surface aabbs
    aabb_t world_box{boxes, boxes.size(), pt_envelope};

    // translation
    const point3 center = world_box.template center<point3>();

    // The world box local frame is the global coordinate frame
    const point3 box_min = world_box.template loc_min<point3>();
    const point3 box_max = world_box.template loc_max<point3>();

    // Get the half lengths for the rectangle sides and translation
    const point3 h_lengths = 0.5f * (box_max - box_min);
    const scalar h_x{math_ns::abs(h_lengths[0])};
    const scalar h_y{math_ns::abs(h_lengths[1])};
    const scalar h_z{math_ns::abs(h_lengths[2])};

    // Volume links for the portal descriptors and the masks
    const dindex volume_idx{volume.index()};
    const nav_link_t volume_link{detail::invalid_value<nav_link_t>()};

    // Construct portals in the...

    //
    // ... x-y plane
    //
    // Only one rectangle needed for both surfaces
    mask_link_t mask_link{rectangle_id, masks.template size<rectangle_id>()};
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_x,
                                              h_y);

    // No rotation, but shift in z for both faces
    vector3 shift{0.f, 0.f, std::isinf(h_z) ? max_shift : h_z};
    transforms.emplace_back(ctx, static_cast<vector3>(center + shift));
    transforms.emplace_back(ctx, static_cast<vector3>(center - shift));

    // Add material slab (no material on the portals)
    material_link_t material_link{slab_id, materials.template size<slab_id>()};
    materials.template emplace_back<slab_id>(empty_context{}, vacuum<scalar>{},
                                             0.f);

    // Build the portal surfaces
    dindex trf_idx{transforms.size(ctx) - 2};
    surfaces.emplace_back(trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);

    surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);

    //
    // ... x-z plane
    //
    ++mask_link;
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_x,
                                              h_z);

    // Rotate by 90deg around x-axis, plus shift in y
    shift = {0.f, std::isinf(h_y) ? max_shift : h_y, 0.f};
    vector3 new_x{1.f, 0.f, 0.f};
    vector3 new_z{0.f, -1.f, 0.f};
    transforms.emplace_back(ctx, static_cast<vector3>(center + shift), new_z,
                            new_x);
    transforms.emplace_back(ctx, static_cast<vector3>(center - shift), new_z,
                            new_x);

    ++material_link;
    materials.template emplace_back<slab_id>(empty_context{}, vacuum<scalar>{},
                                             0.f);

    surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);

    surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);

    //
    // ... y-z plane
    //
    ++mask_link;
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_z,
                                              h_y);

    // Rotate by 90deg around y-axis, plus shift in x
    shift = {std::isinf(h_x) ? max_shift : h_x, 0.f, 0.f};
    new_x = {0.f, 0.f, -1.f};
    new_z = {1.f, 0.f, 0.f};
    transforms.emplace_back(ctx, static_cast<vector3>(center + shift), new_z,
                            new_x);
    transforms.emplace_back(ctx, static_cast<vector3>(center - shift), new_z,
                            new_x);

    ++material_link;
    materials.template emplace_back<slab_id>(empty_context{}, vacuum<scalar>{},
                                             0.f);

    surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);

    surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                          dindex_invalid, surface_id::e_portal);
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
template <auto mask_id, typename context_t, typename trajectory_t,
          typename volume_t, typename surface_container_t,
          typename mask_container_t, typename material_container_t,
          typename transform_container_t, typename config_t>
inline void create_telescope(context_t &ctx, const trajectory_t &traj,
                             volume_t &volume, surface_container_t &surfaces,
                             mask_container_t &masks,
                             material_container_t &materials,
                             transform_container_t &transforms,
                             const config_t &cfg) {
    // @Note: These type definitions should be removed at some point
    using vector3 = __plugin::vector3<detray::scalar>;

    using surface_type = typename surface_container_t::value_type;
    using nav_link_t = typename surface_type::navigation_link;
    using mask_link_t = typename surface_type::mask_link;
    using material_link_t = typename surface_type::material_link;

    auto volume_idx = volume.index();
    constexpr auto slab_id = material_link_t::id_type::e_slab;

    // Create the module centers
    const std::vector<module_placement> m_placements =
        module_positions(traj, cfg.dists);

    // Create geometry data
    for (const auto &m_placement : m_placements) {

        auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

        // Surfaces with the linking into the local containers
        mask_link_t mask_link{mask_id, masks.template size<mask_id>()};
        material_link_t material_link{slab_id,
                                      materials.template size<slab_id>()};
        const auto trf_index = transforms.size(ctx);
        surfaces.emplace_back(trf_index, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_sensitive);

        // The rectangle bounds for this module
        masks.template emplace_back<mask_id>(empty_context{}, cfg.mask_values,
                                             mask_volume_link);

        // Lines need different material
        using mask_t = typename mask_container_t::template get_type<mask_id>;
        if (mask_t::shape::name == "line") {
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
    covfie::field<
        typename telescope_types<typename mask_t::shape>::bfield_backend_t>
        &&bfield,
    const mask_t &msk, std::vector<scalar> dists,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um,
    const trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    const scalar envelope = 0.1f * unit<scalar>::mm) {

    // detector type
    using detector_t = detector<telescope_types<typename mask_t::shape>,
                                covfie::field, container_t>;

    // @todo: Temporal restriction due to missing local navigation
    assert((dists.size() < 20u) &&
           "Due to WIP, please choose less than 20 surfaces for now");

    // module parameters
    struct surface_config {
        const typename mask_t::mask_values &mask_values;
        const std::vector<scalar> &dists;
        material<scalar> m_mat;
        scalar m_thickness;
    };

    // create empty detector
    detector_t det(resource, std::move(bfield));

    typename detector_t::geometry_context ctx{};
    const surface_config sf_config{msk.values(), dists, mat, thickness};

    // Dummy volume bounds for now, will be set correctly when portals are built
    auto &vol = det.new_volume(volume_id::e_cuboid);

    // Add module surfaces to volume
    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    constexpr auto mask_id{
        detector_t::mask_container::template get_id<mask_t>::value};

    create_telescope<mask_id>(ctx, traj, vol, surfaces, masks, materials,
                              transforms, sf_config);
    // Add portals to volume
    create_cuboid_portals<mask_id>(ctx, envelope, vol, surfaces, masks,
                                   materials, transforms);
    // Add volme to the detector
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
          typename container_t = host_container_types>
auto create_telescope_detector(
    vecmem::memory_resource &resource,
    covfie::field<
        typename telescope_types<typename mask_t::shape>::bfield_backend_t>
        &&bfield,
    const mask_t &msk, std::size_t n_surfaces = 10u,
    const scalar tel_length = 500.f * unit<scalar>::mm,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um,
    const trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    const scalar envelope = 0.01f * unit<scalar>::mm) {
    // Generate equidistant positions
    std::vector<scalar> distances = {};
    scalar pos = 0.f;
    scalar dist{n_surfaces > 1u
                    ? tel_length / static_cast<scalar>(n_surfaces - 1u)
                    : 0.f};
    for (std::size_t i = 0u; i < n_surfaces; ++i) {
        distances.push_back(pos);
        pos += dist;
    }

    // Build the geometry
    return create_telescope_detector<mask_t>(resource, std::move(bfield), msk,
                                             distances, mat, thickness, traj,
                                             envelope);
}

/// Wrapper for create_telescope_geometry with constant zero bfield.
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types>
auto create_telescope_detector(
    vecmem::memory_resource &resource, const mask_t &msk,
    const std::vector<scalar> dists,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um,
    const trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    const scalar envelope = 0.1f * unit<scalar>::mm) {

    using covfie_bkdn_t =
        typename telescope_types<typename mask_t::shape>::bfield_backend_t;

    // Build the geometry
    return create_telescope_detector<mask_t, trajectory_t, container_t>(
        resource,
        covfie::field<covfie_bkdn_t>{
            typename covfie_bkdn_t::configuration_t{0.f, 0.f, 0.f}},
        msk, dists, mat, thickness, traj, envelope);
}

/// Wrapper for create_telescope_geometry with constant zero bfield.
template <typename mask_t = mask<rectangle2D<>>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>,
          typename container_t = host_container_types>
auto create_telescope_detector(
    vecmem::memory_resource &resource, const mask_t &msk,
    std::size_t n_surfaces = 10u,
    const scalar tel_length = 500.f * unit<scalar>::mm,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80.f * unit<scalar>::um,
    const trajectory_t traj = {{0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, -1.f},
    const scalar envelope = 0.1f * unit<scalar>::mm) {

    using covfie_bkdn_t =
        typename telescope_types<typename mask_t::shape>::bfield_backend_t;

    // Build the geometry
    return create_telescope_detector<mask_t, trajectory_t, container_t>(
        resource,
        covfie::field<covfie_bkdn_t>{
            typename covfie_bkdn_t::configuration_t{0.f, 0.f, 0.f}},
        msk, n_surfaces, tel_length, mat, thickness, traj, envelope);
}

}  // namespace detray
