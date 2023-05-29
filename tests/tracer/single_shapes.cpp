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
#include "detray/masks/sphere2D.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracer/renderer/camera.hpp"
#include "detray/tracer/renderer/detail/mask.hpp"
#include "detray/tracer/renderer/pipeline.hpp"
#include "detray/tracer/renderer/raw_image.hpp"
#include "detray/tracer/shaders/background.hpp"
#include "detray/tracer/shaders/material.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

// System include(s)
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>

using namespace detray;

namespace {

/// Simple color gradient
template <typename color_depth>
inline void write_test_image(raw_image<color_depth> &im) {
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
template <typename color_depth, typename aspect_ratio, typename mask_t,
          typename material_t, class im_background_t = gradient_background>
inline void render_single_shape(raw_image<color_depth, aspect_ratio> &im,
                                const mask_t &mask, const transform3D &trf,
                                const material_t &mat) {
    using scalar_t = transform3D::scalar_type;

    // Rendering steps
    using intersector_t = single_shape<mask_t, material_t>;
    using backgr_shader_t = background_shader<inf_plane<im_background_t>>;
    using mat_shader_t = material_shader;
    // The rendering pipeline: The intersector finds the shape intersections
    using pipeline_t =
        rendering_pipeline<intersector_t, backgr_shader_t, mat_shader_t>;

    const scalar_t viewport_height{2.0f};
    const point3D origin{0.0f, 0.0f, 0.0f};

    camera<scalar_t, aspect_ratio> cam(viewport_height, origin);

    // For the single shape render, the scene is actually encoded directly in
    // the single shape intersector
    typename intersector_t::global_state geo{trf, mask, mat};

    // Iterate through pixel matrix
    for (std::size_t i_y = 0u; i_y < im.height(); ++i_y) {
        for (std::size_t i_x = 0u; i_x < im.width(); ++i_x) {

            // Ray to render the pixel at (i_x, i_y)
            detail::ray<transform3D> ray = cam.get_ray(i_x, i_y, im);
            ray.set_overstep_tolerance(-std::numeric_limits<scalar>::max());

            // Strap the global geometry state and the thread-local ray together
            scene_handle::state scene{geo, im, ray, i_x, i_y};

            // Finds the intersections between the ray and the geometry
            typename intersector_t::state intrs{};
            actor::state empty{};  // Dummy for actors that don't define a state
            auto pipeline_state = std::tie(empty, intrs);

            // Run
            pipeline_t{}(pipeline_state, scene);

            im.set_pixel(scene.m_pixel);
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
    raw_image<> image{500u};
    write_test_image(image);
    ppm.write(image, "test");

    //
    // Render single shape
    //

    // Affine transform matrix to place the shapes
    vector3D x{1.0f, 0.0f, 0.0f};
    vector3D z{0.0f, 0.0f, 1.f};
    vector3D t{0.0f, 0.0f, 30.0f};
    transform3D trf{t, z, x};

    const silicon_tml<scalar> sf_mat{};

    // render a rectangle mask
    /*const mask<rectangle2D<>> rect{0u, 12.f, 20.f};
    render_single_shape(image, rect, trf, beryllium<scalar>{});
    ppm.write(image, "rectangle");

    // render a trapezoid mask
    const mask<trapezoid2D<>> trpz{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    render_single_shape<>(image, trpz, trf, aluminium<scalar>{});
    ppm.write(image, "trapezoid");

    // render a ring mask
    const mask<ring2D<>> ring{0u, 12.f, 20.f};
    render_single_shape<>(image, ring, trf, gold<scalar>{});
    ppm.write(image, "ring");

    // render an annulus mask
    const mask<annulus2D<>> ann2{0u,       5.f,  13.0f, 0.74195f,
                                 1.33970f, -2.f, 2.f,   0.f};
    render_single_shape<>(image, ann2, trf, silicon<scalar>{});
    ppm.write(image, "annulus");*/

    // render a spherical mask
    const tracer_mask<sphere2D<>> sph2{0u, 10.f};
    render_single_shape<>(image, sph2, trf, silicon<scalar>{});
    ppm.write(image, "sphere");

    return EXIT_SUCCESS;
}
