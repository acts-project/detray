
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {
/** A zone finder based on a grid
 *
 * @tparam grid_t is the type of the underlying grid object
 */
template <typename grid_t>
struct local_zone_finder {
    grid_t _grid;
    bool _sort = true;

    using point2 = __plugin::point2<detray::scalar>;

    /** Constructor from grid
     *
     * @param grid is the prepared object/index grid, it will be moved into the
     * the local finder object
     **/
    local_zone_finder(grid_t &&grid) : _grid(std::move(grid)) {}

    /** Call operator for the object search with binned neighborhood
     *
     * @param p2 the local 2d point for the grid
     * @param nhood the local binned neighborhood
     *
     * This method will create a vector and fill it
     *
     * @note return a binned zone around a bin
     **/
    auto operator()(const point2 &p2,
                    const typename grid_t::template neighborhood<dindex>
                        &nhood = grid_t::hermit2) const {
        return _grid.zone(p2, nhood, _sort);
    }

    /** Call operator for the object search with scalar neighborhood
     *
     * @param p2 the local 2d point for the grid
     * @param nhood the local scalar neighborhood
     *
     * This method will create a vector and fill it
     *
     * @note return a binned zone around a bin
     **/
    auto operator()(
        const point2 &p2,
        const typename grid_t::template neighborhood<scalar> &nhood) const {
        return _grid.zone(p2, nhood, _sort);
    }

    /** Const access to the grid */
    const grid_t &grid() const { return _grid; }

    /** Non-const access to the grid */
    grid_t &grid() { return _grid; }
};

/** A zone finder for a single object */
template <typename value_t, typename vector_t = dvector<value_t>,
          template <typename, std::size_t> class array_t = darray>
struct local_single_finder {

    using point2 = __plugin::point2<detray::scalar>;

    vector_t _value = {};

    /** Constructor from a single value */
    local_single_finder(value_t &&value) { _value = {std::move(value)}; }

    /** Call operator for the object search
     *
     * @tparam point2_type the type of the point for the finding request
     * @param p2 the local 2d point for the grid
     * @note return always the same bin
     **/
    vector_t operator()(const point2 & /*p2*/,
                        const array_t<unsigned int, 2> & /*nhood*/ = {
                            0, 0}) const {
        return _value;
    }
};

}  // namespace detray