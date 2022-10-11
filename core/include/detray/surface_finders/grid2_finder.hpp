/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"

namespace detray {

template <template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple>
using regular_circular_grid =
    grid2<attach_populator, axis::regular, axis::circular, serializer2,
          vector_t, jagged_vector_t, array_t, tuple_t, dindex, false>;

template <template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple>
using regular_circular_grid2 =
    grid2<replace_populator, axis::regular, axis::circular, serializer2,
          vector_t, jagged_vector_t, array_t, tuple_t, dindex, false>;

template <template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple>
using regular_regular_grid =
    grid2<attach_populator, axis::regular, axis::regular, serializer2, vector_t,
          jagged_vector_t, array_t, tuple_t, dindex, false>;

}  // namespace detray