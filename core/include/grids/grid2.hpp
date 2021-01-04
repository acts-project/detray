/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/containers.hpp"

namespace detray
{

    /** A two-dimensional grid for object storage
     * 
     * @tparam populator_type  is a prescription what to do when a bin gets pupulated, it broadcasts
     *         also the value type
     * @tparam tparam axis_p0_type the type of the first axis
     * @tparam tparam axis_p1_type the type of the second axis
     * @tparam serialzier_type  type of the serializer to the storage represenations
     * 
     **/
    template <typename populator_type,
              typename axis_p0_type,
              typename axis_p1_type,
              typename serializer_type>
    class grid2
    {

        axis_p0_type _axis_p0;
        axis_p1_type _axis_p1;
        serializer_type _serializer;
        populator_type _populator;
        darray<typename populator_type::store_value, axis_p0_type::bins * axis_p1_type::bins> _data_serialized;

        /** Constructor from axes (moved)
         * 
         * @param axis_p0 is the axis in the first coordinate
         * @param axis_p1 is the axis in the second coordinate
         * 
         **/
        grid2(axis_p0_type &&ap0, axis_p1_type &&ap1) : _axis_p0(std::move(_axis_p0)), _axis_p1(std::move(_axis_p1)) {}

        /** Fill/populate operation
         * 
         * @param fvalue is a single fill value to be filled
         * @param p2 the point in p2 local frame
         * 
         **/
        template <typename point2_type>
        void populate(typename populator_type::bare_value &&fvalue, const point2_type &p2)
        {
            auto sbin = _serializer.serialize({_axis_p0.bin(p2[0]), _axis_p1.bin(p2[1])});
            _data_serialized[sbin] = _populator(_data_serialized[sbin], std::move(fvalue));
        }

        /** Return the value of a single bin 
         * 
         * @param p2 is the local coordinate p2 type 
         * 
         * @return the const reference to the value in this bin 
         **/
        template <typename point2_type>
        const auto &bin(const point2_type &p2) const
        {
            return _data_serialized[_serializer.serialize({_axis_p0.bin(p2[0]), _axis_p1.bin(p2[1])})];
        }

        /** Return the value of a single bin - non-const access
         * 
         * @param p2 is the local coordinate p2 type 
         * 
         * @return the const reference to the value in this bin 
         **/
        template <typename point2_type>
        auto &bin(const point2_type &p2)
        {
            return _data_serialized[_serializer.serialize({_axis_p0.bin(p2[0]), _axis_p1.bin(p2[1])})];
        }

        template <typename point2_type>
        dvector<typename populator_type::bare_value> zone(const point2_type &p2, const darray<unsigned int, 2> &nhood = {0, 0}) const
        {
            return {};
        }
    };

} // namespace detray