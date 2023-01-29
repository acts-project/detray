/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_array_definitions.hpp"

// Project include(s).
#include "detray/intersection/intersection.hpp"
#include "detray/io/image/ppm_writer.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracer/color.hpp"

// System include(s)
#include <iostream>

using namespace detray;

int main() {

    constexpr color<> red{139u, 0u, 0u, 0u};
    constexpr color<> blue{0u, 0u, 139u, 0u};
    constexpr color<> purple{red + blue};

    io::raw_image out_im{100u, 100u, purple};

    io::ppm_writer ppm("test");
    ppm.write(out_im);

    return 0;
}
