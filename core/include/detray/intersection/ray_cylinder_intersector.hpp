/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/quadratic_equation.hpp"

namespace detray {

namespace detail {

struct unbound;

}

/// @brief Intersection implementation for cylinder surfaces using ray
/// trajectories.
struct ray_cylinder_intersector {

    using intersection_type = line_plane_intersection;

    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using cylindrical2 = __plugin::cylindrical2<detray::scalar>;

    /// Intersection method for cylindrical surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam track_t The type of the track carrying the context object
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param track the track information
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename track_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) const {

        return intersect(trf, track.pos(), track.dir(), mask, tolerance,
                         track.overstep_tolerance());
    }

    /// Intersection method for cylindrical surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param ro the origin of the ray
    /// @param rd the direction of the ray
    /// @param local to the local local frame
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    /// @param overstep_tolerance  is the stepping specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const point3 &ro, const vector3 &rd,
        const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon,
        const scalar overstep_tolerance = 0.) const {

        using local_frame = typename mask_t::local_type;

        const scalar r{mask[0]};
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0, 2);
        const vector3 sc = getter::vector<3>(m, 0, 3);

        const vector3 pc_cross_sz = vector::cross(ro - sc, sz);
        const vector3 rd_cross_sz = vector::cross(rd, sz);
        const scalar a = vector::dot(rd_cross_sz, rd_cross_sz);
        const scalar b = scalar{2.} * vector::dot(rd_cross_sz, pc_cross_sz);
        const scalar c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

        quadratic_equation<scalar> qe = {a, b, c};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > scalar{0.}) {
            const auto t01 = std::get<1>(qe_solution);
            const scalar t{(t01[0] > overstep_tolerance) ? t01[0] : t01[1]};
            if (t > overstep_tolerance) {
                intersection_type is;
                is.path = t;
                is.p3 = ro + is.path * rd;
                constexpr local_frame local_converter{};
                is.p2 = local_converter(trf, is.p3);
                const auto local3 = trf.point_to_local(is.p3);
                is.status =
                    mask.template is_inside<local_frame>(local3, tolerance);
                const scalar rdr = getter::perp(local3 + scalar{0.1} * rd);
                is.direction = rdr > r ? intersection::direction::e_along
                                       : intersection::direction::e_opposite;
                is.link = mask.volume_link();
                return is;
            }
        }
        return intersection_type{};
    }
};

}  // namespace detray
