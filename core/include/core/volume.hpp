/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

namespace detray
{
    /** Templated volume class
     * 
     * @tparam surface_container_type type of the surface container
     * @tparam portal_container_type the type of the portal container 
     * @tparam source_type the type of the volume representation 
     */
    template <typename surface_container_type,
              typename portal_container_type = int,
              typename source_type = int>
    struct volume
    {
        surface_container_type _surfaces;
        portal_container_type _portals;
        source_type _source;

        /** Constructor from arguments
         * 
         * @param surfaces the container of surfaces in this volume, includind boundary surfaces
         * @param portals the container of portals to other volumes
         * @param source the source object described by this volume
         **/
        volume(surface_container_type &&surfaces, portal_container_type &&portals, source_type &&_source)
            : _surfaces(std::move(surfaces)),
              _portals(std::move(portals)),
              _source(std::move(source))
        {
        }

        /** Retrieve the surface candidates from this volume
         * 
         * @param p point at which the query takes place
         * @param d direction of the queary track/particle
         * @param b_field the magnetic field vector 
         * @param mom the absolute momentum at the query
         * 
         * @note boundary surfaces are excluded in the return type
         * 
         * @return an iterable container
         **/
        template <typename point_type, typename vector_type>
        auto surface_candidates(const point_type &p, const vector_type &d, const vector_type &b_field, scalar mom = 0.) const
        {
            return _surfaces.get_candidates(p, d, b_field, mom);
        }

        /** Retrieve the portal candidates from this volume
         * 
         * @param p point at which the query takes place
         * @param d direction of the queary track/particle
         * @param b_field the magnetic field vector 
         * @param mom the absolute momentum at the query
         * 
         * @return an iterable container
         **/
        template <typename point_type, typename vector_type>
        auto portal_candidates(const point_type &p, const vector_type &d, const vector_type &b_field, scalar mom = 0.) const
        {
            return _portals.get_candidates(p, d, b_field, mom);
        }
    };
} // namespace detray