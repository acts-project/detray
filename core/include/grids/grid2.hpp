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

    template <typename value_type, 
              typename axis_p0_type, 
              typename axis_p1_type, 
              typename serializer_type>
    struct grid2
    {
        darray<value_type, typename axis_p0_type::kBins, * typename axis_p1_type::kBins> _data_serialized;


        value_type bin() const {

        }

        dvector<value_type> zone() const {


        }

    };

} // namespace detray