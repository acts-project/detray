/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/accessor.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/tools/grid_array_helper.hpp"

#pragma once

namespace detray {

/** The surfaces finder definition.
 *
 * This class is a surfaces finder definition which wraps around the array of
 * grid that plays a role of local navigation
 *
 * @tparam N the number of grid elements
 * @tparam array_type the type of the internal array, must have STL semantics
 * @tparam tuple_type the type of the internal tuple, must have STL semantics
 * @tparam vector_type the type of the vector, must have STL semantics
 * @tparam jagged_vector_type the type of the vector of vector, must have STL
 * semantics
 */
template <std::size_t N,
          template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          template <typename...> class jagged_vector_type = djagged_vector>
struct surfaces_finder {

    static constexpr size_t N_GRIDS = N;

    template <typename T, unsigned int I>
    using array_t = array_type<T, I>;

    template <typename... Args>
    using tuple_t = tuple_type<Args...>;

    using surfaces_regular_circular_grid =
        grid2<attach_populator, axis::regular, axis::circular, serializer2,
              vector_type, jagged_vector_type, array_type, tuple_type, dindex,
              false>;

    /** Host costructor
     * @param resource memory resource for the allocation of members
     */
    DETRAY_HOST surfaces_finder(vecmem::memory_resource& resource)
        : _surface_grids(initialize_host_grids(
              resource, std::make_index_sequence<N_GRIDS>{})) {}

    /** Aggregate initialization of host grids
     * @tparam ints index_sequence for aggregate initialization
     * @param resource memory resource for the allocation of members
     * @param seq STL index_sequence object for template parameter deduction
     */
    template <std::size_t... ints>
    DETRAY_HOST auto initialize_host_grids(
        vecmem::memory_resource& resource,
        std::index_sequence<ints...> /*seq*/) {

        // TODO: Are raw pointers OK here?
        array_type<vecmem::memory_resource*, N_GRIDS> resources;

        std::fill(resources.begin(), resources.end(), &resource);
        return array_type<surfaces_regular_circular_grid, N_GRIDS>(
            {{*resources[ints]...}});
    }

    /** Device costructor
     * @tparam surfaces_finder_data_type The type of input surface_finder_data
     * object
     * @param finder_data input surface_finder_data object
     */
    template <typename surfaces_finder_data_type,
              std::enable_if_t<!std::is_base_of_v<vecmem::memory_resource,
                                                  surfaces_finder_data_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE surfaces_finder(surfaces_finder_data_type& finder_data)
        : _surface_grids(initialize_device_grids(
              finder_data, std::make_index_sequence<N_GRIDS>{})) {}

    /** Aggregate initialization of device grids
     * @tparam ints index_sequence for aggregate initialization
     * @tparam surfaces_finder_data_type The type of input surface_finder_data
     * object
     * @param finder_data input surface_finder_data object
     * @param seq STL index_sequence object for template parameter deduction
     */
    template <std::size_t... ints, typename surfaces_finder_data_type>
    DETRAY_HOST_DEVICE auto initialize_device_grids(
        surfaces_finder_data_type& finder_data,
        std::index_sequence<ints...> /*seq*/) {
        return array_type<surfaces_regular_circular_grid, N_GRIDS>(
            {{finder_data._surface_grids_view[ints]...}});
    }

    /** return the size of array of grids
     * @return the size of array of grids
     */
    DETRAY_HOST_DEVICE constexpr size_t size() const { return N_GRIDS; }

    /** Access operator - non-const
     * @return the surface grid element
     */
    DETRAY_HOST_DEVICE
    auto& operator[](unsigned int value_index) {
        return _surface_grids[value_index];
    }

    /** Array of grids for the local surface navigation  **/
    array_type<surfaces_regular_circular_grid, N_GRIDS> _surface_grids;
};

/** A static implemetation of surfaces_finder data for device
 * @tparam surfaces_finder_type The type of host surface_finder
 */
template <typename surfaces_finder_type>
struct surfaces_finder_data {
    // type definition
    using surface_grid_t =
        typename surfaces_finder_type::surfaces_regular_circular_grid;

    surfaces_finder_data(surfaces_finder_type& finder,
                         vecmem::memory_resource& resource)
        : _surface_grids_data(
              make_grid_data_array(finder._surface_grids, resource)) {}

    typename surfaces_finder_type::template array_t<
        grid2_data<surface_grid_t>, surfaces_finder_type::N_GRIDS>
        _surface_grids_data;
};

/** A static implemetation of surfaces_finder view for device
 * @tparam surfaces_finder_type The type of host surface_finder
 */
template <typename surfaces_finder_type>
struct surfaces_finder_view {
    // type definition
    using surface_grid_t =
        typename surfaces_finder_type::surfaces_regular_circular_grid;

    surfaces_finder_view(surfaces_finder_data<surfaces_finder_type>& data)
        : _surface_grids_view(make_grid_view_array(data._surface_grids_data)) {}

    typename surfaces_finder_type::template array_t<
        grid2_view<surface_grid_t>, surfaces_finder_type::N_GRIDS>
        _surface_grids_view;
};

/** standalone function for surface_finder_data get function
 */
template <size_t N, template <typename, unsigned int> class array_type,
          template <typename...> class tuple_type,
          template <typename...> class vector_type,
          template <typename...> class jagged_vector_type>
inline surfaces_finder_data<
    surfaces_finder<N, array_type, tuple_type, vector_type, jagged_vector_type>>
get_data(surfaces_finder<N, array_type, tuple_type, vector_type,
                         jagged_vector_type>& finder,
         vecmem::memory_resource& resource) {
    return {finder, resource};
}

}  // namespace detray