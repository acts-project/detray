/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <detray/grids/grid2.hpp>
#include <vecmem/containers/static_array.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

template <typename grid_type, std::size_t N, std::size_t... ints>
auto make_grid_data_array(vecmem::static_array<grid_type, N>& arr,
                          vecmem::memory_resource& resource,
                          std::index_sequence<ints...> seq) {

    return vecmem::static_array<grid2_data<grid_type>, N>{
        get_data(arr[ints], resource)...};
}

template <typename grid_type, std::size_t N>
auto make_grid_data_array(vecmem::static_array<grid_type, N>& arr,
                          vecmem::memory_resource& resource) {
    auto seq = std::make_index_sequence<N>{};
    return make_grid_data_array(arr, resource, seq);
}

template <typename grid_type, std::size_t N, std::size_t... ints>
auto make_grid_view_array(vecmem::static_array<grid2_data<grid_type>, N>& arr,
                          std::index_sequence<ints...> seq) {
    return vecmem::static_array<grid2_view<grid_type>, N>{arr[ints]...};
}

template <typename grid_type, std::size_t N>
auto make_grid_view_array(vecmem::static_array<grid2_data<grid_type>, N>& arr) {
    auto seq = std::make_index_sequence<N>{};
    return make_grid_view_array(arr, seq);
}

}  // namespace detray