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

    public:
        using serialized_storage = dvector<typename populator_type::store_value>;

        /** Constructor from axes (moved)
         * 
         * @param axis_p0 is the axis in the first coordinate
         * @param axis_p1 is the axis in the second coordinate
         * 
         **/
        grid2(axis_p0_type &&axis_p0, axis_p1_type &&axis_p1) : _axis_p0(std::move(axis_p0)), _axis_p1(std::move(axis_p1))
        {
            _data_serialized = serialized_storage(axis_p0.bins() * axis_p1.bins(), _populator.init());
        }

        /** Allow for grid shift, when using a centralized store and indices
        * 
        * @param offset is the applied offset shift
        * 
        **/
        void shift(const typename populator_type::bare_value &offset)
        {
            std::for_each(_data_serialized.begin(), _data_serialized.end(), [&](auto &ds) { _populator.shift(ds, offset); });
        }

        /** Fill/populate operation
         * 
         * @tparam point2_type the 2D local point type
         * 
         * @param p2 the point in p2 local frame
         * @param fvalue is a single fill value to be filled
         * 
         **/
        template <typename point2_type>
        void populate(const point2_type &p2, typename populator_type::bare_value &&fvalue)
        {
            auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]));
            _populator(_data_serialized[sbin], std::move(fvalue));
        }

        /** Return the value of a single bin 
         * 
         * @tparam point2_type the 2D local point type
         * 
         * @param p2 is point in the local frame
         * 
         * @return the const reference to the value in this bin 
         **/
        template <typename point2_type>
        const auto &bin(const point2_type &p2) const
        {
            return _data_serialized[_serializer.template serialize<axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
        }

        /** Return the value of a single bin - non-const access
         * 
         * @tparam point2_type the 2D local point type
         * 
         * @param p2 is point in the local frame
         * 
         * @return the const reference to the value in this bin 
         **/
        template <typename point2_type>
        auto &bin(const point2_type &p2)
        {
            return _data_serialized[_serializer.template serialize<axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
        }

        /** Return a zone around a single bin 
         * 
         * The zone is done with a symmetric neighborhood around the bin which is defined by p2
         *          
         * @tparam point2_type the 2D local point type
         * 
         * @param p2 is point in the local frame
         * @param nhood is the bin-wise neighborhood
         * @param sort is a directive whether to sort or not
         * 
         * @return the sequence of values
         **/
        template <typename point2_type>
        dvector<typename populator_type::bare_value> zone(const point2_type &p2, const darray<unsigned int, 2> &nhood = {0, 0}, bool sort = false) const
        {
            auto zone0 = _axis_p0.zone(p2[0], nhood[0]);
            auto zone1 = _axis_p1.zone(p2[1], nhood[1]);

            dvector<typename populator_type::bare_value> zone;

            // Specialization for bare value equal to store value
            if constexpr (std::is_same_v<typename populator_type::bare_value, typename populator_type::store_value>)
            {
                unsigned int iz = 0;
                zone = dvector<typename populator_type::bare_value>(zone0.size() * zone1.size(), {});
                for (const auto z1 : zone1)
                {
                    for (const auto z0 : zone0)
                    {
                        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, z0, z1);
                        zone[iz++] =  _data_serialized[sbin];
                    }
                }
            }
            else
            {
                zone.reserve(10);
                for (const auto z1 : zone1)
                {
                    for (const auto z0 : zone0)
                    {
                        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, z0, z1);
                        auto bin_data = _data_serialized[sbin];
                        auto bin_content = _populator.sequence(bin_data);
                        zone.insert(zone.end(), bin_content.begin(), bin_content.end());
                    }
                }
            }

            if (sort)
            {
                std::sort(zone.begin(), zone.end());
            }
            return zone;
        }

        /** Const access to axis p0  */
        const axis_p0_type &axis_p0() const { return _axis_p0; }

        /** Const access to axis p1 */
        const axis_p1_type &axis_p1() const { return _axis_p1; }

        /* Copy of axes in a tuple */
        dtuple<axis_p0_type, axis_p1_type> axes() const { return std::tie(_axis_p0, _axis_p1); }

        /** Const acess to the serializer */
        const serializer_type &serializer() const { return _serializer; }

        /** Const acess to the polulator */
        const populator_type &populator() const { return _populator; }

    private:
        serialized_storage _data_serialized;
        axis_p0_type _axis_p0;
        axis_p1_type _axis_p1;
        serializer_type _serializer;
        populator_type _populator;
    };

} // namespace detray