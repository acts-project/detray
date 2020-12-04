/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <vector>
#include <map>

namespace detray
{
    template <typename value_type, unsigned int kDIM>
    using darray = std::array<value_type, kDIM>;

    template <typename value_type>
    using dvector = std::vector<value_type>;

    template <typename key_type, typename value_type>
    using dmap = std::map<key_type, value_type>;

} // namespace detray
