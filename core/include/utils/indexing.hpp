/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

namespace detray
{    
    using dindex = unsigned long;
    auto dindex_invalid = std::numeric_limits<dindex>::max();
    using dindex_range = darray<dindex, 2>;
    using dindex_sequence = dvector<dindex>;

} // namespace detray