/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_array_definitions.hpp"

// Project include(s).
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/io/image/ppm_writer.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/predefined_materials.hpp"
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
inline void write_test_image(io::raw_image<color_depth> &im) {
    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0.f}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
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
template <typename color_depth, typename shape_t, typename material_t,
          template <typename> class im_background_t = gradient_background>
inline void render_single_shape(io::raw_image<color_depth> &im,
                                const detray::mask<shape_t> &mask,
                                const transform3D &trf, const material_t &mat) {
    // Rendering steps
    using intersector_t = single_shape<detray::mask<shape_t>, material_t>;
    using colorizer_t = colorizer<im_background_t>;
    // The full pipeline
    using renderer_t = composite_actor<dtuple, intersector_t, colorizer_t>;
    // @todo: add shaders that observe the renderer ...
    using pipeline_t = rendering_pipeline<renderer_t>;

    const point3D lower_left_corner{-10.0f, -5.0f, -5.0f};
    const point3D horizontal{20.0f, 0.0f, 0.0f};
    const point3D vertical{0.0f, 10.0f, 0.0f};
    const point3D origin{0.0f, 0.0f, 0.0f};

    // For the single shape render, the scene is actually encoded directly in
    // the single shape intersector
    typename intersector_t::global_state geo{trf, mask, mat};

    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
            const float u{i_x / im.width()};
            const float v{i_y / im.height()};

            // Ray to render the pixel at (i_x, i_y)
            detail::ray<transform3D> ray{
                origin, 0.f, lower_left_corner + u * horizontal + v * vertical,
                0.f};

            // Strap the global geometry state and the thread-local ray together
            const scene_handle::state scene{geo, ray};

            // Finds the intersections between the ray and the geometry
            typename intersector_t::state intrs{};
            // Computes an initial color for the ray's pixel
            typename colorizer_t::state colr{static_cast<std::uint64_t>(i_x), static_cast<std::uint64_t>(i_y)};
            // Feed their thread-local state to the pipeline
            auto pipeline_state = std::tie(intrs, colr);

            // Run
            pipeline_t{}(pipeline_state, scene);

            im.set_pixel(colr.m_pixel);
        }
    }
}

}  // namespace

int main() {

    io::ppm_writer<> ppm{};

    //
    // Test image
    //

    // write a test image
    io::raw_image<> image{1000u, 500u};
    write_test_image(image);
    ppm.write(image, "test");

    //
    // Render single shape
    //

    // Affine transform matrix to place the shapes
    vector3D x{1.0f, 0.0f, 0.0f};
    vector3D z{0.0f, 0.0f, 1.f};
    vector3D t{5.0f, 5.0f, 30.0f};
    const transform3D trf{t, z, x};

    const silicon_tml<scalar> sf_mat{};

    // render a rectangle mask
    const mask<rectangle2D<>> rect{0u, 12.f, 20.f};
    render_single_shape(image, rect, trf, sf_mat);
    ppm.write(image, "rectangle");

    // render a trapezoid mask
    const mask<trapezoid2D<>> trpz{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    render_single_shape<>(image, trpz, trf, sf_mat);
    ppm.write(image, "trapezoid");

    // render a ring mask
    const mask<ring2D<>> ring{0u, 12.f, 20.f};
    render_single_shape<>(image, ring, trf, sf_mat);
    ppm.write(image, "ring");

    // render an annulus mask
    const mask<annulus2D<>> ann2{0u,       5.f,  13.0f, 0.74195f,
                                 1.33970f, -2.f, 2.f,   0.f};
    render_single_shape<>(image, ann2, trf, sf_mat);
    ppm.write(image, "annulus");

    return EXIT_SUCCESS;
}
