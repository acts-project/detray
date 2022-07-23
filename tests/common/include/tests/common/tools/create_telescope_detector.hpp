/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "tests/common/tools/detector_metadata.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <limits>

namespace detray {

namespace {

using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

using telescope_types = detector_registry::telescope_detector;

struct module_placement {
    point3 _pos;
    vector3 _dir;

    bool operator==(const module_placement &other) const {
        bool is_same = true;
        is_same &= (_pos == other._pos);
        is_same &= (_dir == other._dir);
        return is_same;
    }
};

/// Helper method for positioning the plane surfaces
///
/// @param track pilot track along which the modules should be placed
/// @param b_field determines the trajectory
/// @param tel_length the length of the telescope
/// @param n_surfaces the number of plane surfaces
///
/// @return a vector of the module positions along the trajectory
template <typename track_t, typename stepper_t>
inline std::vector<module_placement> module_positions(
    track_t &track, stepper_t &stepper, std::vector<scalar> &steps) {
    // dummy navigation struct
    struct navigation_state {
        scalar operator()() const { return _step_size; }
        inline void set_full_trust() {}
        inline void set_high_trust() {}
        inline void set_fair_trust() {}
        inline void set_no_trust() {}
        inline bool abort() { return false; }

        scalar _step_size = 0;
    };

    // dummy propagator state
    struct prop_state {
        typename stepper_t::state _stepping;
        navigation_state _navigation;
    };

    // create and fill the positions
    std::vector<module_placement> m_positions;
    m_positions.reserve(steps.size());

    // Find exact position by walking along track
    prop_state propagation{typename stepper_t::state{track},
                           navigation_state{}};

    // Calculate step size from module positions. The modules will only be
    // placed at the given position if the b-field allows for it. Otherwise, by
    // virtue of the stepper, the distance will be shorter and the surface
    // remains reachable in a single step.
    scalar prev_dist = 0.;
    for (const auto &dist : steps) {
        // advance the track state to the next plane position
        propagation._navigation._step_size = dist - prev_dist;
        stepper.step(propagation);
        m_positions.push_back(
            {propagation._stepping().pos(), propagation._stepping().dir()});
        prev_dist = dist;
    }

    return m_positions;
}

/// Helper function that creates the telescope surfaces.
///
/// @param mask_id id of the plane surface shape, either unbounded or recangular
/// @param ctx geometric context
/// @param track pilot track along which the modules should be placed
/// @param b_field determines the trajectory
/// @param volume volume the planes should be added to
/// @param surfaces container to add new surface to
/// @param masks container to add new cylinder mask to
/// @param transforms container to add new transform to
/// @param cfg config struct for module creation
template <telescope_types::mask_ids mask_id, typename context_t,
          typename track_t, typename stepper_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename material_container_t, typename transform_container_t,
          typename config_t>
inline void create_telescope(context_t &ctx, track_t &track, stepper_t &stepper,
                             volume_type &volume, surface_container_t &surfaces,
                             mask_container_t &masks,
                             material_container_t &materials,
                             transform_container_t &transforms, config_t &cfg) {
    using surface_type = typename surface_container_t::value_type;
    using volume_link_t = typename surface_type::volume_link_type;
    using mask_link_type = typename surface_type::mask_link;
    using material_defs = typename surface_type::material_defs;
    using material_link_type = typename surface_type::material_link;

    auto volume_id = volume.index();
    volume_link_t mask_volume_link{volume_id};
    constexpr auto slab_id = material_defs::id::e_slab;

    // Create the module centers
    const std::vector<module_placement> m_placements =
        module_positions(track, stepper, cfg.pos);

    // Create geometry data
    for (const auto &m_placement : m_placements) {

        // Surfaces with the linking into the local containers
        mask_link_type mask_link{mask_id, masks.template size<mask_id>()};
        material_link_type material_link{slab_id,
                                         materials.template size<slab_id>()};
        const auto trf_index = transforms.size(ctx);
        surfaces.emplace_back(trf_index, mask_link, material_link, volume_id,
                              dindex_invalid, false);
        surfaces.back().set_grid_status(false);

        // The last surface acts as portal that leaves the telescope
        if (m_placement == m_placements.back()) {
            mask_volume_link = dindex_invalid;
        }

        if constexpr (mask_id ==
                      telescope_types::mask_ids::e_unbounded_plane2) {
            // No bounds for this module
            masks.template add_value<
                telescope_types::mask_ids::e_unbounded_plane2>(
                mask_volume_link);
            materials.template add_value<telescope_types::material_ids::e_slab>(
                cfg.m_mat, cfg.m_thickness);
        } else {
            // The rectangle bounds for this module
            masks.template add_value<telescope_types::mask_ids::e_rectangle2>(
                cfg.m_half_x, cfg.m_half_y, mask_volume_link);
            materials.template add_value<telescope_types::material_ids::e_slab>(
                cfg.m_mat, cfg.m_thickness);
        }
        // Build the transform
        // Local z axis is the global normal vector
        vector3 m_local_z = algebra::vector::normalize(m_placement._dir);

        // Project onto the weakest direction component of the normal vector
        vector3 e_i{0., 0., 0.};
        scalar min = std::numeric_limits<scalar>::infinity();
        dindex i = dindex_invalid;
        for (unsigned int k = 0; k < 3; ++k) {
            if (m_local_z[k] < min) {
                min = m_local_z[k];
                i = k;
            }
        }
        e_i[i] = 1.;
        vector3 proj = algebra::vector::dot(m_local_z, e_i) * m_local_z;
        // Local x axis is the normal to local y,z
        vector3 m_local_x = algebra::vector::normalize(e_i - proj);

        // Create the global-to-local transform of the module
        transforms.emplace_back(ctx, m_placement._pos, m_local_z, m_local_x);
    }
}

}  // namespace

/// Build the telescope geometry from a fixed length and number of surfaces.
///
/// @param n_surfaces the number of surfaces that are placed in the geometry
/// @param tel_length the total length of the steps by the stepper
template <bool unbounded_planes = true,
          typename track_t = free_track_parameters,
          typename stepper_t = line_stepper<track_t>,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector>
auto create_telescope_detector(
    vecmem::memory_resource &resource, dindex n_surfaces = 10,
    scalar tel_length = 500. * unit_constants::mm,
    track_t track = {{0., 0., 0.}, 0., {0., 0., 1.}, -1.},
    stepper_t stepper = {}, scalar half_x = 20. * unit_constants::mm,
    scalar half_y = 20. * unit_constants::mm) {
    // Generate equidistant positions
    std::vector<scalar> positions = {};
    scalar pos = 0.;
    scalar dist = tel_length / (n_surfaces - 1);
    for (std::size_t i = 0; i < n_surfaces; ++i) {
        positions.push_back(pos);
        pos += dist;
    }

    // Build the geometry
    return create_telescope_detector<unbounded_planes>(
        resource, positions, track, stepper, half_x, half_y);
}

/// Builds a detray geometry that contains only one volume with plane surfaces,
/// where the last surface is the portal that leaves the telescope. The
/// detector is auto-constructed by following a track state through space with
/// the help of a dedicated stepper. The track and stepper determine the
/// positions of the plane surfaces, so that the stepping distance between
/// them is evenly spaced along the length of the telescope.
///
/// @tparam unbounded_planes build the telescope with unbounded plane surfaces
///        (true) or rectangle surfaces (false)
/// @tparam track_t the type of the pilot track
/// @tparam stepper_t the stepper type that advances the track state
///
/// @param resource the memory resource for the detector containers
/// @param pos the module positions. These only correspond to the actual module
///            positions for a straight line pilot track
/// @param track the pilot track along which the surfaces are positioned
/// @param stepper that advances the track through the length of the telescope
/// @param half_x the x half length of the recangle mask
/// @param half_y the y half length of the recangle mask
///
/// @returns a complete detector object
template <bool unbounded_planes = true,
          typename track_t = free_track_parameters,
          typename stepper_t = line_stepper<track_t>,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector>
auto create_telescope_detector(
    vecmem::memory_resource &resource, std::vector<scalar> pos = {},
    track_t track = {{0, 0, 0}, 0, {0, 0, 1}, -1}, stepper_t stepper = {},
    scalar half_x = 20. * unit_constants::mm,
    scalar half_y = 20. * unit_constants::mm,
    const material<scalar> mat = silicon_tml<scalar>(),
    const scalar thickness = 80 * unit_constants::um) {

    // detector type
    using detector_t =
        detector<telescope_types, array_t, tuple_t, vector_t, jagged_vector_t>;

    // module parameters
    struct plane_config {
        scalar m_half_x;
        scalar m_half_y;
        std::vector<scalar> &pos;
        material<scalar> m_mat;
        scalar m_thickness;
    };

    // create empty detector
    detector_t det(resource);

    typename detector_t::context ctx{};
    plane_config pl_config{half_x, half_y, pos, mat, thickness};

    // volume boundaries are not needed. Same goes for portals
    det.new_volume({0., 0., 0., 0., -M_PI, M_PI});
    typename detector_t::volume_type &vol = det.volume_by_index(0);

    // Add module surfaces to volume
    typename detector_t::surface_container surfaces(&resource);
    typename detector_t::mask_container masks = {resource};
    typename detector_t::material_container materials = {resource};
    typename detector_t::transform_container transforms = {resource};

    if constexpr (unbounded_planes) {
        create_telescope<telescope_types::mask_ids::e_unbounded_plane2>(
            ctx, track, stepper, vol, surfaces, masks, materials, transforms,
            pl_config);
    } else {
        create_telescope<telescope_types::mask_ids::e_rectangle2>(
            ctx, track, stepper, vol, surfaces, masks, materials, transforms,
            pl_config);
    }

    det.add_objects_per_volume(ctx, vol, surfaces, masks, materials,
                               transforms);

    return det;
}

}  // namespace detray