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

/// Render a shape
template <typename shape_t, typename color_depth>
inline void render_single_shape(io::raw_image<color_depth>& im,
                                const mask<shape_t>& m) {

    // using intersector_t = single_shape<mask<rectangle2D<>>>;
    using intersector_t = single_shape<mask<shape_t>>;
    using colorizer_t = colorizer<>;
    using renderer_t = composite_actor<dtuple, intersector_t, colorizer_t>;
    // The full pipeline
    using pipeline_t = rendering_pipeline<dtuple, renderer_t>;

    vector3D x{1.0f, 0.0f, 0.0f};
    vector3D z{0.0f, 0.0f, 1.f};
    vector3D t{5.0f, 5.0f, 30.0f};
    transform3D trf{t, z, x};

    const point3D lower_left_corner{-10.0f, -5.0f, -5.0f};
    const point3D horizontal{20.0f, 0.0f, 0.0f};
    const point3D vertical{0.0f, 10.0f, 0.0f};
    const point3D origin{0.0f, 0.0f, 0.0f};

    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
            const float u{i_x / im.width()};
            const float v{i_y / im.height()};

            typename intersector_t::state geo(
                m, trf, origin,
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

    io::ppm_writer<> ppm{};

    // write a test image
    io::raw_image<> image{1000u, 500u};
    write_test_image(image);
    ppm.write(image, "test");

    // render a rectangle mask
    mask<rectangle2D<>> rect{0u, 12.f, 20.f};
    render_single_shape(image, rect);
    ppm.write(image, "rectangle");

    // render a trapezoid mask
    mask<trapezoid2D<>> trpz{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    render_single_shape(image, trpz);
    ppm.write(image, "trapezoid");

    // render a ring mask
    mask<ring2D<>> ring{0u, 12.f, 20.f};
    render_single_shape(image, ring);
    ppm.write(image, "ring");

    // render an annulus mask
    mask<annulus2D<>> ann2{0u, 5.f, 13.0f, 0.74195f, 1.33970f, -2.f, 2.f, 0.f};
    render_single_shape(image, ann2);
    ppm.write(image, "annulus");

    return EXIT_SUCCESS;
}
