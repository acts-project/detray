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
         * The axis is closed, i.e. each underflow bin is mapped to 0
         * and henceforth each overflow bin is mapped to bins-1
         */
        struct regular
        {
            dindex n_bins;
            scalar min;
            scalar max;

            static constexpr unsigned int axis_identifier = 0;

            /** Return the number of bins */
            dindex bins() const { return n_bins; }

            /** Access function to a single bin from a value v
             * 
             * @param v is the value for the bin search
             * 
             * As the axis is closed it @returns a dindex type
             **/
            dindex bin(scalar v) const
            {
                int ibin = static_cast<int>((v - min) / (max - min) * n_bins);
                return (ibin >= 0 and ibin < n_bins) ? static_cast<dindex>(ibin)
                       : ibin < 0                    ? 0
                                                     : static_cast<dindex>(n_bins - 1);
            }

            /** Access function to a range with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a dindex_range
             **/
            dindex_range range(scalar v, unsigned int nhood) const
            {

                int ibin = static_cast<int>((v - min) / (max - min) * n_bins);
                int ibinmin = ibin - static_cast<int>(nhood);
                int ibinmax = ibin + static_cast<int>(nhood);
                dindex min_bin = (ibinmin >= 0) ? static_cast<dindex>(ibinmin) : 0;
                dindex max_bin = (ibinmax < static_cast<int>(n_bins)) ? static_cast<dindex>(ibinmax) : static_cast<dindex>(n_bins - 1);
                return {min_bin, max_bin};
            }

            /** Access function to a zone with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a dindex_sequence
             **/
            dindex_sequence zone(scalar v, unsigned int nhood) const
            {
                dindex_range nh_range = range(v, nhood);
                dindex_sequence sequence(static_cast<dindex_sequence::size_type>(nh_range[1] - nh_range[0] + 1), nh_range[0]);
                dindex m = 0;
                std::for_each(sequence.begin(), sequence.end(), [&](auto &n) { n += m++; });
                return sequence;
            }

            /** Copy the range zone */
            darray<scalar, 2> range() const { return {min, max}; }
        };

        /** A regular circular axis.
         * 
         * The axis is circular, i.e. the underflow bins map into the circular sequence
         */
        struct circular
        {

            unsigned int n_bins;
            scalar min;
            scalar max;

            static constexpr unsigned int axis_identifier = 1;

            /** Return the number of bins */
            dindex bins() const { return n_bins; }

            /** Access function to a single bin from a value v
             * 
             * @param v is the value for the bin search
             * 
             * As the axis is closed it @returns a dindex type
             **/
            dindex bin(scalar v) const
            {
                dindex ibin = static_cast<dindex>((v - min) / (max - min) * n_bins);
                return (ibin >= 0 and ibin < n_bins) ? static_cast<dindex>(ibin)
                       : ibin < 0                    ? static_cast<dindex>(n_bins + ibin)
                                                     : static_cast<dindex>(n_bins - ibin);
            }

            /** Access function to a range with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is circular it @returns a dindex_range
             **/
            dindex_range range(scalar v, unsigned int nhood) const
            {
                dindex gbin = bin(v);
                dindex min_bin = remap(gbin, -static_cast<int>(nhood));
                dindex max_bin = remap(gbin, static_cast<int>(nhood));
                return {min_bin, max_bin};
            }

            /** Access function to a zone with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a dindex_sequence
             **/
            dindex_sequence zone(scalar v, unsigned int nhood) const
            {
                dindex_range nh_range = range(v, nhood);
                if (nh_range[0] < nh_range[1])
                {
                    dindex_sequence sequence(static_cast<dindex_sequence::size_type>(nh_range[1] - nh_range[0] + 1), nh_range[0]);
                    dindex m = 0;
                    std::for_each(sequence.begin(), sequence.end(), [&](auto &n) { n += m++; });
                    return sequence;
                }
                dindex vl = static_cast<dindex>(n_bins - nh_range[0] + nh_range[1] + 1);
                dindex mi = 0;
                dindex mo = 0;
                dindex_sequence sequence(static_cast<dindex_sequence::size_type>(vl), nh_range[0]);
                std::for_each(sequence.begin(), sequence.end(), [&](auto &n) {
                    n += mi++;
                    if (n > n_bins - 1)
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
             * @return an index, remapped bin 
             **/
            dindex remap(dindex ibin, int shood) const
            {
                int opt_bin = static_cast<int>(ibin) + shood;
                if (opt_bin >= 0 and opt_bin < n_bins)
                {
                    return static_cast<dindex>(opt_bin);
                }
                if (opt_bin < 0)
                {
                    return static_cast<dindex>(n_bins + opt_bin);
                }
                return static_cast<dindex>(opt_bin - n_bins);
            }

            /** Copy the range zone */
            darray<scalar, 2> range() const { return {min, max}; }
        };

        /** An iregular circular axis.
         * 
         * The axis is closed, i.e. the underflow is mapped into the first,
         * the overflow is mapped into the last.
         */
        struct irregular
        {

            dvector<scalar> boundaries;

            static constexpr unsigned int axis_identifier = 2;

            /** Return the number of bins */
            dindex bins() const { return static_cast<dindex>(boundaries.size() - 1); }

            /** Access function to a single bin from a value v
             * 
             * @param v is the value for the bin search
             * 
             * As the axis is closed it @returns a dindex type
             **/
            dindex bin(scalar v) const
            {
                return (std::lower_bound(boundaries.begin(), boundaries.end(), v) - boundaries.begin());
            }

            /** Access function to a range with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a dindex_range
             **/
            dindex_range range(scalar v, unsigned int nhood) const
            {

                unsigned int ibin = bin(v);
                int bins = boundaries.size() - 1;
                int ibinmin = ibin - static_cast<int>(nhood);
                int ibinmax = ibin + static_cast<int>(nhood);
                dindex min_bin = (ibinmin >= 0) ? static_cast<dindex>(ibinmin) : 0;
                dindex max_bin = (ibinmax < static_cast<int>(bins)) ? static_cast<dindex>(ibinmax) : static_cast<dindex>(bins - 1);
                return {min_bin, max_bin};
            }

            /** Access function to a zone with binned neighbourhood
             * 
             * @param v is the value for the bin search
             * @param nhood is the neighbourhood size (+/-) 
             * 
             * As the axis is closed it @returns a dindex_sequence
             **/
            dindex_sequence zone(scalar v, unsigned int nhood) const
            {
                dindex_range nh_range = range(v, nhood);
                dindex_sequence sequence(static_cast<dindex_sequence::size_type>(nh_range[1] - nh_range[0] + 1), nh_range[0]);
                dindex m = 0;
                std::for_each(sequence.begin(), sequence.end(), [&](auto &n) { n += m++; });
                return sequence;
            }

            /** Copy the range zone */
            darray<scalar, 2> range() const { return {boundaries[0], boundaries[boundaries.size() - 1]}; }
        };

    } // namespace axis

} // namespace detray