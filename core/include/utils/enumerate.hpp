/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "utils/indexing.hpp"

#include <tuple>

namespace detray
{

    /** Helper utility to allow indexed enumeration with structured binding
     * 
     * Usage:
     * 
     * for (auto [ i, value ] = enumerate(container) ) { ... };
     * 
     * with 'container' any stl-like container
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

    /** Helper method to (fake a) run over a single entry
     * 
     * Usage:
     * for (auto i : sequence(j)) {}
     * 
     * with j an unsinged dindex type
     * 
     * @note sequence(2) will produce {2}
     * 
     **/
    constexpr auto sequence(dindex iterable)
    {

        struct iterator
        {
            size_t i;

            bool operator!=(const iterator &rhs) const { return i != rhs.i; }

            /** Increase index and iterator at once */
            void operator++()
            {
                ++i;
            }

            /** Tie them together for returning */
            auto operator*() const { return i; }
        };
        struct iterable_wrapper
        {
            dindex iterable;
            auto begin() { return iterator{iterable}; }
            auto end() { return iterator{iterable + 1}; }
        };
        return iterable_wrapper{std::forward<dindex>(iterable)};
    }

    /** Helper method to run over a range 
     * 
     * Usage:
     * for (auto i : sequence(r)) {}
     * 
     * with r an range that can be accessed with r[0] and r[1]
     * 
     * @note sequence({2,4}) will produce { 2, 3, 4 }
     * 
     **/
    constexpr auto sequence(darray<dindex, 2> iterable)
    {

        struct iterator
        {
            size_t i;
            size_t end;

            bool operator!=(const iterator &rhs) const { return i != rhs.end; }

            /** Increase index and iterator at once */
            void operator++()
            {
                ++i;
            }

            /** Tie them together for returning */
            auto operator*() const { return i; }
        };
        struct iterable_wrapper
        {
            darray<dindex, 2> iterable;
            auto begin() { return iterator{iterable[0], iterable[1]}; }
            auto end() { return iterator{iterable[1] + 1, iterable[1] + 1}; }
        };
        return iterable_wrapper{std::forward<darray<dindex, 2> >(iterable)};
    }

} // namespace detray
