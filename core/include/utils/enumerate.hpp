/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Licenced under: Apache-2, see LICENSE file
 */
#pragma once

#include <tuple>

namespace detray
{

    /** Helper utility to allow indexed enumeration with structured binding
     * Usage:
     * 
     * auto [ i, value ] = enumerate(container);
     * 
     */
    template <typename container_type,
              typename container_type_iter = decltype(std::begin(std::declval<container_type>())),
              typename = decltype(std::end(std::declval<container_type>()))>
    constexpr auto enumerate(container_type &&iterable)
    {
        struct iterator
        {
            size_t i;
            container_type_iter iter;
            
            bool operator!=(const iterator &rhs) const { return iter != rhs.iter; }

            /** Increase index and iterator at once */
            void operator++()
            {
                ++i;
                ++iter;
            }

            /** Tie them together for returning */
            auto operator*() const { return std::tie(i, *iter); }
        };
        struct iterable_wrapper
        {
            container_type iterable;
            auto begin() { return iterator{0, std::begin(iterable)}; }
            auto end() { return iterator{0, std::end(iterable)}; }
        };
        return iterable_wrapper{std::forward<container_type>(iterable)};
    }

} // namespace detray
