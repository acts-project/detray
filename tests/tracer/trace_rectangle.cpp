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
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

// System include(s)
#include <cstdlib>
#include <iostream>

using namespace detray;

int main() {

    constexpr texture::color<> red{139u, 0u, 0u, 0u};
    constexpr texture::color<> blue{0u, 0u, 139u, 0u};
    constexpr texture::color<> purple{red + blue};

    constexpr texture::pixel<> px{{0u, 0u}, purple};

    std::cout << px << std::endl;

    io::raw_image<> out_im{100u, 100u, purple};

    io::ppm_writer ppm("test");
    ppm.write(out_im);

    return EXIT_SUCCESS;
}
