/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/indexing.hpp"
#include "utils/containers.hpp"

#include <algorithm>

namespace detray
{
    /** Struct that helps taking a two-dimensional binning into
     * a serial binning for data storage.
     * 
     * Serializers allow to create a memory local environment if
     * advantegeous. 
     * 
     **/
    struct serializer2
    {
        /** Create a serial bin from two individual bins 
         * 
         * @tparam faxis_type is the type of the first axis
         * @tparam saxis_type is the type of the second axis
         * 
         * @return a guaranteed_index for the memory storage
         */
        template <typename faxis_type, typename saxis_type>
        guaranteed_index serialize(guaranteed_index fbin, guaranteed_index sbin) const
        {
            guaranteed_index offset = sbin * faxis_type::bins;
            return offset + fbin;
        }

        /** Create a bin tuple from a serialized bin
         * 
         * @tparam faxis_type is the type of the first axis
         * @tparam saxis_type is the type of the second axis
         * 
         * @return a 2-dim array of guaranteed_index 
         */
        template <typename faxis_type, typename saxis_type>
        darray<guaranteed_index, 2> deserialize(guaranteed_index serialbin) const
        {
            guaranteed_index sbin = static_cast<guaranteed_index>(serialbin/faxis_type::bins);
            guaranteed_index fbin = serialbin - sbin * faxis_type::bins;
            return { fbin, sbin };
        }
    };

} // namespace detray
