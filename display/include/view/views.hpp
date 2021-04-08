/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/containers.hpp"

#include <cmath>

namespace detray
{
    using point3 = __plugin::transform3::point3;
    using transform3 = __plugin::transform3;

    using contour = darray<dvector<scalar>, 2>;

    // Single view per module/surface
    struct single_view
    {

        /// A single view operator, no transform is done
        ///
        /// @param vertices the vertices of the surface
        /// @param tf the (ignored) surface transform
        ///
        /// return a 2D contour array
        contour
        operator()(const dvector<point3> &vertices, const transform3& /*tf*/) const
        {

            dvector<scalar> x;
            x.reserve(vertices.size());
            dvector<scalar> y;
            y.reserve(vertices.size());

            for (const auto &v : vertices)
            {
                x.push_back(v[0]);
                y.push_back(v[1]);
            }
            return {x, y};
        }
    };

    /// x-y projection view
    struct global_xy_view
    {
        /// A global xy view operator, no transform is done
        ///
        /// @param vertices the vertices of the surface
        /// @param tf the surface transform 
        ///
        /// @return a 2D contour array
        contour
        operator()(const dvector<point3> &vertices, const transform3 &tf) const
        {
            dvector<scalar> x;
            x.reserve(vertices.size());
            dvector<scalar> y;
            y.reserve(vertices.size());

            for (const auto &v : vertices)
            {
                auto vg = tf.point_to_global(v);
                x.push_back(vg[0]);
                y.push_back(vg[1]);
            }
            return {x, y};
        }
    };


    /// r-z projection view
    struct global_rz_view
    {
        /// A global rz view operator, no transform is done
        ///
        /// @param vertices the vertices of the surface
        /// @param tf the surface transform 
        ///
        /// @return a 2D contour array
        contour
        operator()(const dvector<point3> &vertices, const transform3 &tf) const
        {
            dvector<scalar> r;
            r.reserve(vertices.size());
            dvector<scalar> z;
            z.reserve(vertices.size());

            for (const auto &v : vertices)
            {
                auto vg = tf.point_to_global(v);
                r.push_back(getter::perp(vg));
                z.push_back(vg[2]);
            }
            return {z, r};
        }
    };

    // rphi-z projection view
    struct rphi_z_view
    {

        scalar fixed_r = std::numeric_limits<scalar>::quiet_NaN();

        /// A view for a rphi-z projection
        ///
        /// @param vertices the vertices of the surface
        /// @param tf the surface transform 
        ///
        /// @return a 2D contour array 
        contour
        operator()(const dvector<point3> &vertices, const transform3 &tf) const
        {
            dvector<scalar> x;
            x.reserve(vertices.size());
            dvector<scalar> y;
            y.reserve(vertices.size());

            for (const auto &v : vertices)
            {
                auto vg = tf.point_to_global(v);
                scalar r = std::isnan(fixed_r) ? getter::perp(vg) : fixed_r;
                x.push_back(r*getter::phi(vg));
                y.push_back(vg[2]);
            }
            return {x, y};
        }
    };

}