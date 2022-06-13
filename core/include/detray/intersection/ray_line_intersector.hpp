/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

/** This is an intersector struct for line surface
 * The intersectino point is obtained by calculating a transform3 object of
 * virtual plane Virtual plane has a normal vector same with the direction of
 * the track
 */
struct ray_line_intersector {

    using intersection_type = line_plane_intersection;
    using transform3 = __plugin::transform3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /** Function to obtain the virtual plane whose normal vector is same with
     * direction of the track
     *
     * @param trf the surface to be intersected
     * @param track the track information
     *
     * @return the transform3 of the virtual plane
     */
    template <typename track_t>
    transform3 line_to_planar_transform3(const transform3 &trf,
                                         const track_t &track) const {

        // z direction of virtual plane parallel to the track direction
        const vector3 new_z = track.dir();

        // y direction of virtual plane identical to line direction
        const vector3 new_y = getter::vector<3>(trf.matrix(), 0, 2);

        // x direction of virtual plane
        const vector3 new_x = vector::cross(new_y, new_z);

        scalar y_dist =
            std::abs(vector::dot(trf.translation() - track.pos(), new_y));

        const vector3 new_t = trf.translation() - y_dist * new_y;

        return transform3{new_t, new_x, new_y, new_z};
    }

    /** Intersection method for line surfaces
     *
     * @tparam track_t The type of the track (which carries the context
     *         object)
     * @tparam mask_t The mask type applied to the local frame
     * @tparam local_frame The local frame type to be intersected
     *
     * Contextual part:
     * @param trf the surface to be intersected
     * @param track the track information
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <typename track_t, typename mask_t,
              std::enable_if_t<std::is_class_v<typename mask_t::local_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform3 &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) const {

        // Convert line transform3 to planar transform3
        const transform3 planar_trf = line_to_planar_transform3(trf, track);

        // return planar intersector result
        return ray_plane_intersector().intersect(planar_trf, track, mask,
                                                 tolerance);
    }
};

}  // namespace detray
