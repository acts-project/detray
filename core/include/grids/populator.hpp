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
    /** A replace populator that swaps whatever current value in the 
     * bin with the new one.
     * 
     * @tparam value_type the type of a single stored object
     * 
     * @note bare_value and store_value are identicial in this case
     **/
    template <typename value_type>
    struct replace_populator
    {

        using bare_value = value_type;
        using store_value = value_type;

        /** Swap the stored value with a new bare value
         * 
         * @param stored the stored value for the population
         * @param bvalue the new value to be added
         **/
        void operator()(store_value &stored, bare_value &&bvalue)
        {
            stored = std::move(bvalue);
        }
    };

    /** A complete populator that adds values to the internal
     * store array until it is completed, ignored afterwards.
     * 
     * @tparam value_type the type of a single stored object
     * 
     * @note bare_value and store_value are different in this case
     **/
    template <typename value_type, value_type kInvalid, unsigned int kDIM, bool kSORT = false>
    struct complete_populator
    {

        static constexpr value_type invalid_value = kInvalid;

        using bare_value = value_type;
        using store_value = darray<bare_value, kDIM>;

        /** Complete the stored value with a new bare value
         * 
         * @param stored the stored value for the population
         * @param bvalue the new value to be added
         **/
        void operator()(store_value &stored, bare_value &&bvalue)
        {

            for (auto &val : stored)
            {
                if (val == kInvalid)
                {
                    val = bvalue;
                    break;
                }
            }
            if (kSORT)
            {
                std::sort(stored.begin(), stored.end());
            }
        }
    };

    /** An attach populator that adds the new value to the 
     * 
     * @tparam value_type the type of a single stored object
     * 
     * @note bare_value and store_value are identicial in this case
     **/
    template <typename value_type, bool kSORT = false>
    struct attach_populator
    {

        using bare_value = value_type;
        using store_value = dvector<bare_value>;

        /** Add a new value to the stored value
         * 
         * @param stored the stored value for the population
         * @param bvalue the new value to be added
         **/
        void operator()(store_value &stored, bare_value &&bvalue)
        {
            stored.push_back(bvalue);
            if (kSORT)
            {
                std::sort(stored.begin(), stored.end());
            }
        }
    };

} // namespace detray