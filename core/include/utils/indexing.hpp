/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/containers.hpp"

namespace detray
{

    using optional_index = int;
    using guaranteed_index = unsigned long;
    using guaranteed_range = darray<guaranteed_index, 2>;
    using guaranteed_sequence = dvector<guaranteed_index>;

} // namespace detray