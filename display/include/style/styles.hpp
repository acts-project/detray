
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

    using color = darray<float, 4>;

    struct style {
        color fill_color = { 0., 1., 0., 0.};
        scalar line_width = 0.;
        unsigned int segments = 72;
    };

}