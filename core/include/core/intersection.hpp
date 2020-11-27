#pragma once

#include <climits>
#include <optional>

namespace detray
{
    /** This templated class holds the intersection information
     * 
     * @tparam scalar_type is the type of the used scalar for intersecting
     * @tparam gvector_type is the type of global intsersection vector
     * @tparam lvector_type is the type of the local intersection vector
     * 
     **/
    template <typename scalar_type, typename gvector_type, typename lvector_type>
    struct intersection
    {
        scalar_type path = std::numeric_limits<scalar_type>::infinity();
        std::optional<gvector_type> point3 = std::nullopt;
        std::optional<lvector_type> point2 = std::nullopt;

        /** @param rhs is the right hand side intersection for comparison 
         **/
        bool operator<(
            const intersection<scalar_type, gvector_type, lvector_type> &rhs) const
        {
            return (path < rhs.path);
        }

        /** @param rhs is the left hand side intersection for comparison 
         **/
        bool operator>(
           const intersection<scalar_type, gvector_type, lvector_type> &rhs) const
        {
            return (path > rhs.path);
        }
    };

} // namespace detray