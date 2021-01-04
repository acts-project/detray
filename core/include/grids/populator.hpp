/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/containers.hpp"

#include <algorithm>

namespace detray
{
    /** A replace populator that swaps whatever current value in the 
     * bin with the new one.
     * 
     * @tparam value_type the type of a single stored object
     * 
     * @note bare_value and store_value are identicial in this case
     **/
    template <typename value_type, value_type kInvalid>
    struct replace_populator
    {

        using bare_value = value_type;
        using store_value = value_type;

        static constexpr value_type invalid_value = kInvalid;

        /** Swap the stored value with a new bare value
         * 
         * @param stored the stored value for the population
         * @param bvalue the new value to be added
         **/
        void operator()(store_value &stored, bare_value &&bvalue) const
        {
            stored = std::move(bvalue);
        }

        /** Create a sequence of bare values, independent of the store_value.
         * 
         * @param stored the stored value
         * 
         * @return a sequence of bare values
         */
        dvector<bare_value> sequence(store_value &stored) const
        {
            if (stored != kInvalid)
            {
                return {stored};
            }
            return {};
        }

        /** Return an initialized bin value
         */
        store_value init() const
        {
            return kInvalid;
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
        void operator()(store_value &stored, bare_value &&bvalue) const
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

        /** Create a sequence of bare values, independent of the store_value.
         * 
         * @param stored the stored value
         * 
         * @return a sequence of bare values, @note it will ignore invalid entries
         */
        dvector<bare_value> sequence(store_value &stored) const
        {
            dvector<bare_value> s;
            s.reserve(kDIM);
            for (const auto &val : stored)
            {
                if (val != kInvalid)
                {
                    s.push_back(val);
                }
            }
            return s;
        }

        /** Return an initialized bin value
         **/
        store_value init() const
        {

            store_value init_bin;
            for (auto &val : init_bin)
            {
                val = kInvalid;
            }
            return init_bin;
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
        void operator()(store_value &stored, bare_value &&bvalue) const
        {
            stored.push_back(bvalue);
            if (kSORT)
            {
                std::sort(stored.begin(), stored.end());
            }
        }

        /** Create a sequence of bare values, independent of the store_value.
         * 
         * @param stored the stored value
         * 
         * @return a sequence of bare values
         */
        dvector<bare_value> sequence(store_value &stored) const
        {
            return stored;
        }

        /** Return an initialized bin value
         **/
        store_value init() const
        {
            return {};
        }
    };

} // namespace detray