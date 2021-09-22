/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>

#include "utils/indexing.hpp"

namespace detray {
/** Struct that helps taking a two-dimensional binning into
 * a serial binning for data storage.
 *
 * Serializers allow to create a memory local environment if
 * advantegeous.
 *
 **/
struct serializer2 {
    /** Create a serial bin from two individual bins
     *
     * @tparam faxis_type is the type of the first axis
     * @tparam saxis_type is the type of the second axis
     *
     * @param faxis the first axis
     * @param saxis the second axis, unused here
     * @param fbin first bin
     * @param sbin second bin
     *
     * @return a dindex for the memory storage
     */
    template <typename faxis_type, typename saxis_type>
    dindex serialize(const faxis_type &faxis, const saxis_type & /*saxis*/,
                     dindex fbin, dindex sbin) const {

        dindex offset = sbin * faxis.bins();
        return offset + fbin;
    }

    /** Create a bin tuple from a serialized bin
     *
     * @tparam faxis_type is the type of the first axis
     * @tparam saxis_type is the type of the second axis
     *
     * @param faxis the first axis
     * @param saxis the second axis, unused here
     * @param serialbin the serial bin
     *
     * @return a 2-dim array of dindex
     */
    template <typename faxis_type, typename saxis_type,
              template <typename, unsigned int> class array_type = darray>
    array_type<dindex, 2> deserialize(const faxis_type &faxis,
                                      const saxis_type & /*saxis*/,
                                      dindex serialbin) const {
        dindex sbin = static_cast<dindex>(serialbin / faxis.bins());
        dindex fbin = serialbin - sbin * faxis.bins();
        return {fbin, sbin};
    }
};

}  // namespace detray
