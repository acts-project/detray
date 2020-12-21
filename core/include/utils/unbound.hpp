
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray
{

    /** A representation that is not bound to a local frame */
    struct unbound
    {
        using point2 = int;

        /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame,
              * including the contextual transform into the local 3D frame
              * 
              * @tparam the type of the surface from which also point3 and context type can be deduced
              * 
              */
        template <typename transform_type>
        const std::optional<point2> operator()(const transform_type & /*ignored*/,
                                               const typename transform_type::point3 & /*ignored*/) const
        {
            return std::nullopt;
        }

        /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             */
        template <typename point3_type>
        const std::optional<point2> operator()(const point3_type & /*ignored*/) const
        {
            return std::nullopt;
        }
    };
} // namespace detray