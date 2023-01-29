/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_array_definitions.hpp"

// Project include(s).
#include "detray/masks/masks.hpp"
#include "detray/tracer/color.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/io/image/ppm_writer.hpp"

// System include(s)
#include <iostream>

using namespace detray;

int main() {

    constexpr color<> red{139.f, 0.f, 0.f, 0.f};

    image out_im{100u, 100u, red};

    ppm_writer ppm("all_red");
    ppm.write(out_im);

    return 0;
}
