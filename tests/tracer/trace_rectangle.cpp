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

namespace {

/// Simple color gradient
template <typename color_depth>
inline void write_test_image(io::raw_image<color_depth>& im) {
    // Iterate through pixel matrix
    for (int i_y{im.height() - 1u}; i_y >= 0; --i_y) {
        for (int i_x{0}; i_x < im.width(); ++i_x) {
            float r{static_cast<float>(i_x) / static_cast<float>(im.width())};
            float g{static_cast<float>(i_y) / static_cast<float>(im.height())};
            float b{0.2f};
            uint8_t i_r{static_cast<uint8_t>(255.99f * r)};
            uint8_t i_g{static_cast<uint8_t>(255.99f * g)};
            uint8_t i_b{static_cast<uint8_t>(255.99f * b)};
            im.set_pixel(i_x, i_y, texture::color<>{i_r, i_g, i_b, 0u});
        }
    }
}

}  // namespace

int main() {

    constexpr texture::color<> red{139u, 0u, 0u, 0u};
    constexpr texture::color<> blue{0u, 0u, 139u, 0u};
    constexpr texture::color<> purple{red + blue};
    ;

    io::raw_image<> out_im{200u, 100u};
    write_test_image(out_im);

    io::ppm_writer ppm("test");
    ppm.write(out_im);

    return EXIT_SUCCESS;
}
