/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <map>
#include <tuple>
#include <vecmem/containers/vector.hpp>
#include <vector>

#include "vecmem/containers/jagged_vector.hpp"

namespace detray {
template <typename value_t, unsigned int kDIM>
using darray = std::array<value_t, kDIM>;

template <typename value_t>
using dvector = vecmem::vector<value_t>;

template <typename value_t>
using djagged_vector = vecmem::jagged_vector<value_t>;

template <typename key_t, typename value_t>
using dmap = std::map<key_t, value_t>;

template <class... types>
using dtuple = std::tuple<types...>;

}  // namespace detray
