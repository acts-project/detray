
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/** A representation that is not bound to a local frame */
struct unbound {
    using transform3 = __plugin::transform3;
    using point3 = __plugin::point3;
    using point2 = __plugin::point2;

    /** This method transform from a point from the global 3D cartesian frame to
     *  the local 2D cartesian frame, including the contextual transform into
     *  the local 3D frame
     **/
    const point2 operator()(const transform3& /*ignored*/,
                            const point3& /*ignored*/) const {
        return {std::numeric_limits<scalar>::infinity(),
                std::numeric_limits<scalar>::infinity()};
    }

    /** This method transform from a point from the global 3D cartesian frame to
     *  the local 2D cartesian frame
     */
    const point2 operator()(const point3& /*ignored*/) const {
        return {std::numeric_limits<scalar>::infinity(),
                std::numeric_limits<scalar>::infinity()};
    }
};
}  // namespace detray