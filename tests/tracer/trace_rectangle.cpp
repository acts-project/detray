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
#include "detray/tracer/renderer/pipeline.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

// System include(s)
#include <cstdlib>
#include <iostream>
#include <tuple>

using namespace detray;

namespace {

/// Simple color gradient
template <typename color_depth>
inline void write_test_image(io::raw_image<color_depth>& im) {
    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
            const float r{i_x / im.width()};
            const float g{i_y / im.height()};
            const float b{0.2f};

            const texture::color<> c_grad{static_cast<uint8_t>(255.99f * r),
                                          static_cast<uint8_t>(255.99f * g),
                                          static_cast<uint8_t>(255.99f * b),
                                          0u};

            im.set_pixel(i_x, i_y, c_grad);
        }
    }
}

/// Render a rectangle
template <typename color_depth>
inline void render_rectangle(io::raw_image<color_depth>& im) {

    using intersector_t = single_shape<mask<rectangle2D<>>>;
    using colorizer_t = colorizer<>;
    using renderer_t = composite_actor<dtuple, intersector_t, colorizer_t>;
    // The full pipeline
    using pipeline_t = rendering_pipeline<dtuple, renderer_t>;

    mask<rectangle2D<>> rect{0u, 5.f, 8.f};

    const point3D lower_left_corner{-2.0f, -1.0f, -1.0f};
    const point3D horizontal{4.0f, 0.0f, 0.0f};
    const point3D vertical{0.0f, 2.0f, 0.0f};
    const point3D origin{0.0f, 0.0f, 0.0f};

    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
            const float u{i_x / im.width()};
            const float v{i_y / im.height()};

            intersector_t::state geo(
                rect, origin,
                lower_left_corner + u * horizontal + v * vertical);
            colorizer_t::state col(i_x, i_y);
            auto pipeline_state = std::tie(geo, col);

            // Propagator state
            struct empty_prop_state {};
            empty_prop_state prop_state{};

            // Run
            pipeline_t run_pipeline{};
            run_pipeline(pipeline_state, prop_state);

            im.set_pixel(col.m_pixel);
        }
    }
}

}  // namespace

int main() {

    constexpr texture::color<> red{139u, 0u, 0u, 0u};
    constexpr texture::color<> blue{0u, 0u, 139u, 0u};
    constexpr texture::color<> purple{red + blue};

    io::raw_image<> out_im{200u, 100u};
    // write_test_image(out_im);
    render_rectangle(out_im);

    io::ppm_writer ppm("test");
    ppm.write(out_im);

    return EXIT_SUCCESS;
}
