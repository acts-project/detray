/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "tests/common/tools/detector_metadata.hpp"

namespace detray {

namespace {

using telescope_types = detector_registry::telescope_detector;

constexpr auto rectangle_id = telescope_types::mask_ids::e_rectangle2;
constexpr auto unbounded_id = telescope_types::mask_ids::e_unbounded_plane2;

/** Helper method for positioning the plane surfaces
 *
 * @param track pilot track along which the modules should be placed
 * @param b_field determines the trajectory
 * @param tel_length the length of the telescope
 * @param n_surfaces the number of plane surfaces
 *
 * @return a vector of the module positions along the trajectory
 */
template <typename track_t, typename stepper_t>
inline std::vector<point3> module_positions(track_t &track, stepper_t &stepper,
                                            scalar tel_length,
                                            dindex n_surfaces) {
    // dummy navigation struct
    struct navigation_state {
        scalar operator()() const { return _step_size; }
        inline void set_full_trust() {}
        inline void set_high_trust() {}
        inline void set_fair_trust() {}
        inline void set_no_trust() {}
        inline bool abort() { return false; }

        scalar _step_size;
    };

    // create and fill the positions
    std::vector<point3> m_positions;
    m_positions.reserve(n_surfaces);

    // space between surfaces
    navigation_state n_state{tel_length / (n_surfaces - 1)};

    // Find exact position by walking along track
    typename stepper_t::state step_state(track);

    m_positions.push_back(track.pos());
    for (size_t i = 1; i < size_t{n_surfaces}; ++i) {
        // advance the track state to the next plane position
        stepper.step(step_state, n_state);
        m_positions.push_back(track.pos());
    }

    return m_positions;
}

/** Helper function that creates the telescope surfaces.
 *
 * @param mask_id id of the plane surface shape, either unbounded or recangular
 * @param ctx geometric context
 * @param track pilot track along which the modules should be placed
 * @param b_field determines the trajectory
 * @param volume volume the planes should be added to
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param cfg config struct for module creation
 */
template <telescope_types::mask_ids mask_id, typename context_t,
          typename track_t, typename stepper_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename transform_container_t, typename config_t>
inline void create_telescope(context_t &ctx, track_t &track, stepper_t &stepper,
                             volume_type &volume, surface_container_t &surfaces,
                             mask_container_t &masks,
                             transform_container_t &transforms, config_t cfg) {
    using surface_t = typename surface_container_t::value_type::value_type;
    using edge_t = typename surface_t::edge_type;
    using mask_link_t = typename mask_container_t::link_type;

    auto volume_id = volume.index();
    edge_t mask_edge{volume_id, dindex_invalid};

    // Create the module centers
    const std::vector<point3> m_positions =
        module_positions(track, stepper, cfg.tel_length, cfg.n_surfaces);

    // Create geometry data
    for (const auto &m_center : m_positions) {

        // Surfaces with the linking into the local containers
        mask_link_t mask_link{mask_id, masks.template size<mask_id>()};
        const auto trf_index = transforms[mask_id].size(ctx);
        surfaces[mask_id].emplace_back(trf_index, mask_link, volume_id,
                                       dindex_invalid, false);
        surfaces[mask_id].back().set_grid_status(false);

        // The last surface acts as portal that leaves the telescope
        if (m_center == m_positions.back()) {
            mask_edge = {dindex_invalid, dindex_invalid};
        }

        if constexpr (mask_id == unbounded_id) {
            // No bounds for this module
            masks.template add_mask<unbounded_id>(mask_edge);
        } else {
            // The rectangle bounds for this module
            masks.template add_mask<rectangle_id>(cfg.m_half_x, cfg.m_half_y,
                                                  mask_edge);
        }
        // Build the transform
        // The local phi
        scalar m_phi = algebra::getter::phi(m_center);
        // Local z axis is the normal vector
        vector3 m_local_z{std::cos(m_phi), std::sin(m_phi), 0.};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi), std::cos(m_phi), 0.};

        // Create the module transform
        transforms[mask_id].emplace_back(ctx, m_center, m_local_z, m_local_x);
    }
}

}  // namespace

/** Builds a detray geometry that contains only one volume with plane surfaces,
 *  where the last surface is the portal that leaves the telescope. The
 *  detector is auto-constructed by following a track state through space with
 *  the help of a dedicated stepper. The track and stepper determine the
 *  positions of the plane surfaces, so that the stepping distance between
 *  them is evenly spaced along the length of the telescope.
 *
 * @tparam unbounded_planes build the telescope with unbounded plane surfaces
 *         (true) or rectangle surfaces (false)
 * @tparam track_t the type of the pilot track
 * @tparam stepper_t the stepper type that advances the track state
 *
 * @param resource the memory resource for the detector containers
 * @param track the pilot track along which the surfaces are positioned
 * @param stepper that advances the track through the length of the telescope
 * @param n_surfaces the number of surfaces that are placed in the geometry
 * @param tel_length the total length of the steps by the stepper
 * @param half_x the x half length of the recangle mask
 * @param half_y the y half length of the recangle mask
 *
 * @returns a complete detector object
 */
template <bool unbounded_planes = true, typename track_t, typename stepper_t,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector>
auto create_telescope_detector(vecmem::memory_resource &resource, track_t track,
                               stepper_t &stepper, dindex n_surfaces = 10,
                               scalar tel_length = 500. * unit_constants::mm,
                               scalar half_x = 20. * unit_constants::mm,
                               scalar half_y = 20. * unit_constants::mm) {

    // detector type
    using detector_t =
        detector<telescope_types, array_t, tuple_t, vector_t, jagged_vector_t>;

    // module parameters
    struct plane_config {
        scalar m_half_x;
        scalar m_half_y;
        scalar tel_length;
        dindex n_surfaces;
    };

    // create empty detector
    detector_t det(resource);

    typename detector_t::context ctx{};
    plane_config pl_config{half_x, half_y, tel_length, n_surfaces};

    // volume boundaries are not needed. Same goes for portals
    det.new_volume(
        {tel_length, tel_length, tel_length, tel_length, -M_PI, M_PI});
    typename detector_t::volume_type &vol = det.volume_by_index(0);

    // Add module surfaces to volume
    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_filling_container transforms = {resource};

    if constexpr (unbounded_planes) {
        create_telescope<unbounded_id>(ctx, track, stepper, vol, surfaces,
                                       masks, transforms, pl_config);
    } else {
        create_telescope<rectangle_id>(ctx, track, stepper, vol, surfaces,
                                       masks, transforms, pl_config);
    }

    det.add_objects(ctx, vol, surfaces, masks, transforms);

    return det;
}

}  // namespace detray