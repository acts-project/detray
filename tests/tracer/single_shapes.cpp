/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_soa_definitions.hpp"

// Project include(s).
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/soa/cylinder_intersector.hpp"
#include "detray/intersection/soa/cylinder_portal_intersector.hpp"
#include "detray/intersection/soa/line_intersector.hpp"
#include "detray/intersection/soa/plane_intersector.hpp"
#include "detray/intersection/soa/sphere_intersector.hpp"
#include "detray/io/image/ppm_writer.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/sphere2D.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracer/renderer/camera.hpp"
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

            const texture::color<color_depth> c_grad{
                static_cast<color_depth>(255.99f * r),
                static_cast<color_depth>(255.99f * g),
                static_cast<color_depth>(255.99f * b), 0u};

            im.set_pixel(i_x, i_y, c_grad);
        }
    }
}

/// Render a shape
template <typename T, template <typename> class algebra_t, typename color_depth,
          typename aspect_ratio, typename mask_t, typename material_t,
          class im_background_t = gradient_background<T, ALGEBRA_PLUGIN>>
inline void render_single_shape(raw_image<color_depth, aspect_ratio> &im,
                                const mask_t &mask,
                                const dtransform3D<algebra_t<T>> &trf,
                                const material_t &mat) {

    // Rendering steps
    using intersector_t = single_shape<T, algebra_t, mask_t, material_t>;
    using backgr_shader_t = background_shader<inf_plane<im_background_t>>;
    using mat_shader_t = material_shader<T, ALGEBRA_PLUGIN>;
    // The rendering pipeline: The intersector finds the shape intersections
    using pipeline_t =
        rendering_pipeline<intersector_t, backgr_shader_t, mat_shader_t>;

    const T viewport_height = 2.0f;
    const dpoint3D<ALGEBRA_PLUGIN<T>> origin{0.0f, 0.0f, 0.0f};

    camera<T, ALGEBRA_PLUGIN, aspect_ratio> cam(viewport_height, origin);

    // For the single shape render, the scene is actually encoded directly in
    // the single shape intersector
    typename intersector_t::global_state geo{trf, mask, mat};

    // Iterate through pixel matrix
    for (std::size_t i_y = 0u; i_y < im.height(); ++i_y) {
        for (std::size_t i_x = 0u; i_x < im.width(); ++i_x) {

            // Ray to render the pixel at (i_x, i_y)
            auto ray = cam.get_ray(i_x, i_y, im);
            ray.set_overstep_tolerance(-std::numeric_limits<T>::max());

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

/// Linear algebra implementation using SoA memory layout
template <typename T>
using algebra_v = detray::vc_soa<T>;

/// Linear algebra implementation using AoS memory layout
template <typename T>
using algebra_s = detray::cmath<T>;

int main() {

    io::ppm_writer<unsigned int> ppm{};

    //
    // Test image
    //

    // write a test image
    raw_image<unsigned int> image{500u};
    write_test_image(image);
    ppm.write(image, "test");

    //
    // Render single shape
    //

    using vector3D_s = dvector3D<algebra_s<scalar>>;
    using vector3D_v = dvector3D<algebra_v<scalar>>;

    // Affine transform matrix to place the shapes
    // AoS
    vector3D_s x_s{1.0f, 0.0f, 0.0f};
    vector3D_s z_s{0.0f, 0.0f, 1.f};
    vector3D_s t_s{0.f, 0.0f, -30.0f};

    dtransform3D<algebra_s<scalar>> trf_s{t_s, z_s, x_s};

    // SoA
    vector3D_v x_v{1.0f, 0.0f, 0.0f};
    vector3D_v z_v{0.0f, 0.0f, 1.f};
    vector3D_v t_v{30.0f, -20.0f, 0.0f};
    t_v[0] = t_v[0].Random();
    t_v[0] = 0.1f * (image.width() * t_v[0] - 0.5f * image.width());
    t_v[1] = t_v[1].Random();
    t_v[1] = 0.1f * (image.height() * t_v[1] - 0.5f * image.height());
    t_v[2] = -30.f;

    dtransform3D<algebra_v<scalar>> trf_v{t_v, z_v, x_v};

    // render a rectangle mask

    // AoS
    const mask<rectangle2D<>> rect2_s{0u, 0.01f * image.width(),
                                      0.01f * image.height()};
    render_single_shape<scalar, algebra_s>(image, rect2_s, trf_s,
                                           beryllium<scalar>{});
    ppm.write(image, "rectangle_AoS");

    // SoA
    const mask<rectangle2D<soa::plane_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        rect2_v{0u, 0.01f * image.width() * dsimd<algebra_v, scalar>{}.Random(),
                0.01f * image.height() * dsimd<algebra_v, scalar>{}.Random()};
    render_single_shape<scalar, algebra_v>(image, rect2_v, trf_v,
                                           beryllium<scalar>{});
    ppm.write(image, "rectangle_SoA");

    // render a trapezoid mask

    // AoS
    const mask<trapezoid2D<>> trpz2_s{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    render_single_shape<scalar, algebra_s>(image, trpz2_s, trf_s,
                                           aluminium<scalar>{});
    ppm.write(image, "trapezoid_AoS");

    // SoA
    const mask<trapezoid2D<soa::plane_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        trap2_v{0u, 0.01f * image.width() * dsimd<algebra_v, scalar>{}.Random(),
                0.03f * image.height() * dsimd<algebra_v, scalar>{}.Random(),
                0.02f * image.height() * dsimd<algebra_v, scalar>{}.Random(),
                1.f / 40.f * 1.f / 1000.f * image.height() *
                    dsimd<algebra_v, scalar>{}.Random()};
    render_single_shape<scalar, algebra_v>(image, trap2_v, trf_v,
                                           aluminium<scalar>{});
    ppm.write(image, "trapezoid_SoA");

    // render a ring mask

    // AoS
    const mask<ring2D<>> ring2_s{0u, 12.f, 20.f};
    render_single_shape<scalar, algebra_s>(image, ring2_s, trf_s,
                                           gold<scalar>{});
    ppm.write(image, "ring_AoS");

    // SoA
    const mask<ring2D<soa::plane_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        ring2_v{0u, 23.f * dsimd<algebra_v, scalar>{}.Random(),
                30.f * dsimd<algebra_v, scalar>{}.Random()};
    render_single_shape<scalar, algebra_v>(image, ring2_v, trf_v,
                                           gold<scalar>{});
    ppm.write(image, "ring_SoA");

    // render an annulus mask

    // AoS
    const mask<annulus2D<>> ann2_s{0u,       5.f,  13.0f, 0.74195f,
                                   1.33970f, -2.f, 2.f,   0.f};
    render_single_shape<scalar, algebra_s>(image, ann2_s, trf_s,
                                           silicon<scalar>{});
    ppm.write(image, "annulus_AoS");

    // SoA
    const auto rand = 2.f * dsimd<algebra_v, scalar>{}.Random();
    const mask<annulus2D<soa::plane_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        ann2_v{0u,
               5.f * rand,
               13.0f * rand,
               0.74195f * rand,
               1.33970f * rand,
               -2.f * rand,
               2.f * rand,
               0.f * rand};
    render_single_shape<scalar, algebra_v>(image, ann2_v, trf_v,
                                           silicon<scalar>{});
    ppm.write(image, "annulus_SoA");

    // render a spherical mask

    // AoS
    const mask<sphere2D<>, std::uint_least16_t, algebra_s<scalar>> sph2_s{0u,
                                                                          10.f};

    render_single_shape<scalar, algebra_s>(image, sph2_s, trf_s,
                                           silicon<scalar>{});
    ppm.write(image, "sphere_AoS");

    // SoA
    const mask<sphere2D<soa::sphere_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        sph2_v{0u, 10.f * dsimd<algebra_v, scalar>{}.Random()};

    render_single_shape<scalar, algebra_v>(image, sph2_v, trf_v,
                                           silicon<scalar>{});
    ppm.write(image, "sphere_SoA");

    // render a line mask

    // AoS
    const mask<line<true>, std::uint_least16_t, algebra_s<scalar>> ln2_s{
        0u, 10.f, std::numeric_limits<scalar>::infinity()};

    render_single_shape<scalar, algebra_s>(image, ln2_s, trf_s, gold<scalar>{});
    ppm.write(image, "line_AoS");

    // SoA
    const mask<line<true, soa::line_intersector>, std::uint_least16_t,
               algebra_v<scalar>>
        ln2_v{0u, 10.f * rand, std::numeric_limits<scalar>::infinity() * rand};

    render_single_shape<scalar, algebra_v>(image, ln2_v, trf_v, gold<scalar>{});
    ppm.write(image, "line_SoA");

    // render a cylinder mask

    // AoS
    const mask<cylinder2D<>, std::uint_least16_t, algebra_s<scalar>> cyl2_s{
        0u, 0.5f * image.height(), 0.5f * image.width(), 0.7f * image.width()};

    render_single_shape<scalar, algebra_s>(image, cyl2_s, trf_s,
                                           silicon<scalar>{});
    ppm.write(image, "cylinder_AoS");

    // SoA
    const mask<cylinder2D<false, soa::cylinder_intersector>,
               std::uint_least16_t, algebra_v<scalar>>
        cyl2_v{0u, 0.5f * image.height() * rand, 0.5f * image.width(),
               0.7f * image.width()};

    render_single_shape<scalar, algebra_v>(image, cyl2_v, trf_v,
                                           silicon<scalar>{});
    ppm.write(image, "cylinder_SoA");

    // render a portal cylinder mask

    // AoS
    const mask<cylinder2D<false, cylinder_portal_intersector>,
               std::uint_least16_t, algebra_s<scalar>>
        pt_cyl2_s{0u, 0.1f * image.height(), 0.5f * image.width(),
                  0.7f * image.width()};

    render_single_shape<scalar, algebra_s>(image, pt_cyl2_s, trf_s,
                                           silicon<scalar>{});
    ppm.write(image, "portal_cylinder_AoS");

    // SoA
    const mask<cylinder2D<false, soa::cylinder_portal_intersector>,
               std::uint_least16_t, algebra_v<scalar>>
        pt_cyl2_v{0u, 0.5f * image.height() * rand, 0.5f * image.width(),
                  0.7f * image.width()};

    render_single_shape<scalar, algebra_v>(image, pt_cyl2_v, trf_v,
                                           silicon<scalar>{});
    ppm.write(image, "portal_cylinder_SoA");

    return EXIT_SUCCESS;
}
