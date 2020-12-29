/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/indexing.hpp"

#include <algorithm>

namespace detray
{

    namespace axis
    {
        /** A regular closed axis.
         * 
         * @tparam kDIM the number of bins
         * @tparam value_type the scalar value which is used
         * 
         * The axis is closed, i.e. each underflow bin is mapped to 0
         * and henceforth each overflow bin is mapped to kDIM-1
         */
        template <unsigned int kDIM, typename value_type = scalar>
        struct closed
        {

            value_type _min;
            value_type _max;

            static constexpr unsigned int bins = kDIM;

            /** Access function to a single bin from a value v
             * 
             * @param v is the value for the bin search
             * 
             * As the axis is closed it @returns a guaranteed_index type
             **/
            guaranteed_index bin(value_type v) const
            {
                optional_index ibin = static_cast<optional_index>((v - _min) / (_max - _min) * kDIM);
                return (ibin >= 0 and ibin < kDIM) ? static_cast<guaranteed_index>(ibin)
                       : ibin < 0                  ? 0
                                                   : static_cast<guaranteed_index>(kDIM - 1);
            }

            /** Access function to a range with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a guaranteed_range
             **/
            guaranteed_range range(value_type v, unsigned int nhood) const
            {
                optional_index ibin = static_cast<optional_index>((v - _min) / (_max - _min) * kDIM);
                guaranteed_index min_bin = (ibin - nhood) >= 0 ? static_cast<guaranteed_index>(ibin - nhood) : 0;
                guaranteed_index max_bin = (ibin + nhood) < kDIM ? static_cast<guaranteed_index>(ibin + nhood) : static_cast<guaranteed_index>(kDIM - 1);
                return {min_bin, max_bin};
            }

            /** Access function to a zone with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a guaranteed_sequence
             **/
            guaranteed_sequence zone(value_type v, unsigned int nhood) const
            {
                guaranteed_range nh_range = range(v, nhood);
                guaranteed_sequence sequence(static_cast<guaranteed_sequence::size_type>(nh_range[1] - nh_range[0] + 1), nh_range[0]);
                guaranteed_index m = 0;
                std::for_each(sequence.begin(), sequence.end(), [&](auto &n) { n += m++; });
                return sequence;
            }
        };

        /** A regular circular axis.
         * 
         * @tparam kDIM the number of bins
         * @tparam value_type the scalar value which is used
         * 
         * The axis is circular, i.e. the underflow bins map into the circular sequence
         * 
         */
        template <unsigned int kDIM, typename value_type = scalar>
        struct circular
        {

            value_type _min;
            value_type _max;

            static constexpr unsigned int bins = kDIM;

            /** Access function to a single bin from a value v
             * 
             * @param v is the value for the bin search
             * 
             * As the axis is closed it @returns a guaranteed_index type
             **/
            guaranteed_index bin(value_type v) const
            {
                optional_index ibin = static_cast<optional_index>((v - _min) / (_max - _min) * kDIM);
                return (ibin >= 0 and ibin < kDIM) ? static_cast<guaranteed_index>(ibin)
                       : ibin < 0                  ? static_cast<guaranteed_index>(kDIM + ibin)
                                                   : static_cast<guaranteed_index>(kDIM - ibin);
            }

            /** Access function to a range with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is circular it @returns a guaranteed_range
             **/
            guaranteed_range range(value_type v, unsigned int nhood) const
            {
                guaranteed_index gbin = bin(v);
                guaranteed_index min_bin = remap(gbin, -nhood);
                guaranteed_index max_bin = remap(gbin, nhood);
                return {min_bin, max_bin};
            }

            /** Access function to a zone with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a guaranteed_sequence
             **/
            guaranteed_sequence zone(value_type v, unsigned int nhood) const
            {
                guaranteed_range nh_range = range(v, nhood);
                if (nh_range[0] < nh_range[1])
                {
                    guaranteed_sequence sequence(static_cast<guaranteed_sequence::size_type>(nh_range[1] - nh_range[0] + 1), nh_range[0]);
                    guaranteed_index m = 0;
                    std::for_each(sequence.begin(), sequence.end(), [&](auto &n) { n += m++; });
                    return sequence;
                }
                guaranteed_index vl = static_cast<guaranteed_index>(kDIM - nh_range[0] + nh_range[1] + 1);
                guaranteed_index mi = 0;
                guaranteed_index mo = 0;
                guaranteed_sequence sequence(static_cast<guaranteed_sequence::size_type>(vl), nh_range[0]);
                std::for_each(sequence.begin(), sequence.end(), [&](auto &n) {
                    n += mi++;
                    if (n > kDIM - 1)
                    {
                        n = mo++;
                    } });
                return sequence;
            }

            /** Helper function to remap onto a circular range 
             * 
             * @param ibin is the optional binning value
             * @param shood is the sided neighbour hood
             * 
             * @return a guaranteed, remapped bin 
             **/
            guaranteed_index remap(guaranteed_index ibin, optional_index shood) const
            {
                optional_index opt_bin = ibin + shood;
                if (opt_bin >= 0 and opt_bin < kDIM - 1)
                {
                    return static_cast<guaranteed_index>(opt_bin);
                }
                if (opt_bin < 0)
                {
                    return static_cast<guaranteed_index>(kDIM + opt_bin);
                }
                return static_cast<guaranteed_index>(opt_bin - kDIM);
            }
        };

    } // namespace axis

} // namespace detray