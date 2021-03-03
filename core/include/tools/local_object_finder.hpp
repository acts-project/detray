
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray
{
    /** A zone finder based on a grid 
     * 
     * @tparam grit_type is the type of the underlying grid object 
     */
    template <typename grid_type>
    struct local_zone_finder
    {
        grid_type _grid;
        bool _sort = true;

        /** Constructor from grid 
         * 
         * @param grid is the prepared object/index grid, it will be moved into the 
         * the local finder object
         **/
        local_zone_finder(grid_type &&grid)
            : _grid(std::move(grid)) {}

        /** Call operator for the object search 
         * 
         * @tparam point2_type the type of the point for the finding request
         * @param p2 the local 2d point for the grid
         * 
         * @note return a zone around a bin
         **/
        template <typename point2_type>
        auto operator()(const point2_type &p2, const darray<unsigned int, 2> &nhood = {0, 0}) const
        {
            return _grid.zone(p2, nhood, _sort);
        }

        /** Const access to the grid */
        const grid_type& grid() const { return _grid; }

        /** Non-const access to the grid */
        grid_type& grid() { return _grid; }

    };

    /** A zone finder for a single object */
    template <typename value_type>
    struct local_single_finder
    {

        dvector<value_type> _value = {};

        /** Constructor from a single value */
        local_single_finder(value_type &&value)
        {
            _value = {std::move(value)};
        }

        /** Call operator for the object search
         * 
         * @tparam point2_type the type of the point for the finding request
         * @param p2 the local 2d point for the grid
         * @note return always the same bin 
         **/
        template <typename point2_type>
        dvector<value_type> operator()(const point2_type &p2, const darray<unsigned int, 2> &nhood = {0, 0}) const
        {
            return _value;
        }

    };

} // namespace detray