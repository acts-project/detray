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
#include <chrono>
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
inline void render_single_shape(
    raw_image<color_depth, aspect_ratio> &im, std::vector<mask_t> &&mask,
    const std::vector<dtransform3D<algebra_t<T>>> &trf,
    const std::vector<material_t> &mat) {

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
    typename intersector_t::global_state geo{std::move(trf), std::move(mask),
                                             std::move(mat)};

// Iterate through pixel matrix
#pragma omp parallel for collapse(2)
    for (std::size_t i_y = 0u; i_y < im.height(); ++i_y) {
        for (std::size_t i_x = 0u; i_x < im.width(); ++i_x) {

            // Ray to render the pixel at (i_x, i_y)
            auto ray = cam.get_ray(i_x, i_y, im);

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
    constexpr std::size_t soa_size{dscalar<algebra_v<scalar>>::size()};

    // Affine transform matrix to place the shapes

    // SoA
    vector3D_v x_v{1.0f, 0.0f, 0.0f};
    vector3D_v z_v{0.0f, 0.0f, 1.f};
    vector3D_v t_v{30.0f, -20.0f, 0.0f};
    t_v[0] = t_v[0].Random();
    t_v[0] = 0.1f * (image.width() * t_v[0] - 0.5f * image.width());
    t_v[1] = t_v[1].Random();
    t_v[1] = 0.1f * (image.height() * t_v[1] - 0.5f * image.height());
    t_v[2] = -30.1f;

    std::vector<dtransform3D<algebra_v<scalar>>> trfs_v;
    trfs_v.emplace_back(t_v, z_v, x_v);

    // AoS
    std::vector<dtransform3D<algebra_s<scalar>>> trfs_s;
    trfs_s.reserve(soa_size);
    for (std::size_t i = 0; i < soa_size; ++i) {
        vector3D_s x_s{x_v[0][i], x_v[1][i], x_v[2][i]};
        vector3D_s z_s{z_v[0][i], z_v[1][i], z_v[2][i]};
        vector3D_s t_s{t_v[0][i], t_v[1][i], t_v[2][i]};

        trfs_s.emplace_back(t_s, z_s, x_s);
    }

    // Different materials per surface
    std::vector<material<scalar>> mat{beryllium<scalar>{}, aluminium<scalar>{},
                                      gold<scalar>{},      silicon<scalar>{},
                                      tungsten<scalar>{},  gold<scalar>{},
                                      aluminium<scalar>{}, silicon<scalar>{}};

    // render a rectangle mask

    // AoS
    mask<rectangle2D<>> rect2_s{0u, 0.01f * image.width(),
                                0.01f * image.height()};
    std::vector<mask<rectangle2D<>>> rect2_vec(soa_size, rect2_s);

    auto start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(rect2_vec), trfs_s,
                                           mat);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "\nRectangle AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "rectangle_AoS");

    // SoA
    std::vector<mask<rectangle2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        rect2_v;
    rect2_v.emplace_back(
        0u, 0.01f * image.width() * dsimd<algebra_v, scalar>{}.Random(),
        0.01f * image.height() * dsimd<algebra_v, scalar>{}.Random());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(rect2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Rectangle SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "rectangle_SoA");

    // render a trapezoid mask

    // AoS
    const mask<trapezoid2D<>> trpz2_s{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    std::vector<mask<trapezoid2D<>>> trpz2_vec(soa_size, trpz2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(trpz2_vec), trfs_s,
                                           mat);

    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nTrapezoid AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "trapezoid_AoS");

    // SoA
    std::vector<mask<trapezoid2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        trap2_v;
    trap2_v.emplace_back(
        0u, 0.01f * image.width() * dsimd<algebra_v, scalar>{}.Random(),
        0.03f * image.height() * dsimd<algebra_v, scalar>{}.Random(),
        0.02f * image.height() * dsimd<algebra_v, scalar>{}.Random(),
        1.f / 40.f * 1.f / 1000.f * image.height() *
            dsimd<algebra_v, scalar>{}.Random());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(trap2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Trapezoid SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "trapezoid_SoA");

    // render a ring mask

    // AoS
    const mask<ring2D<>> ring2_s{0u, 12.f, 20.f};
    std::vector<mask<ring2D<>>> ring2_vec(soa_size, ring2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(ring2_vec), trfs_s,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nRing AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "ring_AoS");

    // SoA
    std::vector<mask<ring2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ring2_v;
    ring2_v.emplace_back(0u, 23.f * dsimd<algebra_v, scalar>{}.Random(),
                         30.f * dsimd<algebra_v, scalar>{}.Random());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(ring2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Ring SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "ring_SoA");

    // render an annulus mask

    // AoS
    const mask<annulus2D<>> ann2_s{0u,       5.f,  13.0f, 0.74195f,
                                   1.33970f, -2.f, 2.f,   0.f};
    std::vector<mask<annulus2D<>>> ann2_vec(soa_size, ann2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(ann2_vec), trfs_s,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nAnnulus AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "annulus_AoS");

    // SoA
    const auto rand = 2.f * dsimd<algebra_v, scalar>{}.Random();
    std::vector<mask<annulus2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ann2_v;
    ann2_v.emplace_back(0u, 5.f * rand, 13.0f * rand, 0.74195f * rand,
                        1.33970f * rand, -2.f * rand, 2.f * rand, 0.f * rand);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(ann2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Annulus SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "annulus_SoA");

    // render a spherical mask

    // AoS
    const mask<sphere2D<>, std::uint_least16_t, algebra_s<scalar>> sph2_s{0u,
                                                                          10.f};
    std::vector<mask<sphere2D<>>> sph2_vec(soa_size, sph2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(sph2_vec), trfs_s,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nSphere AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "sphere_AoS");

    // SoA
    std::vector<mask<sphere2D<soa::sphere_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        sph2_v;
    sph2_v.emplace_back(0u, 10.f * dsimd<algebra_v, scalar>{}.Random());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(sph2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Sphere SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "sphere_SoA");

    // render a line mask

    // AoS
    const mask<line<true>, std::uint_least16_t, algebra_s<scalar>> ln2_s{
        0u, 10.f, std::numeric_limits<scalar>::infinity()};
    std::vector<mask<line<true>>> ln2_vec(soa_size, ln2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(ln2_vec), trfs_s,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nLine AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "line_AoS");

    // SoA
    std::vector<mask<line<true, soa::line_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ln2_v;
    ln2_v.emplace_back(0u, 10.f * rand,
                       std::numeric_limits<scalar>::infinity() * rand);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(ln2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Line SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "line_SoA");

    // render a cylinder mask

    // AoS
    const mask<cylinder2D<>, std::uint_least16_t, algebra_s<scalar>> cyl2_s{
        0u, 0.5f * image.height(), 0.5f * image.width(), 0.7f * image.width()};
    std::vector<mask<cylinder2D<>>> cyl2_vec(soa_size, cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(cyl2_vec), trfs_s,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nCylinder AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "cylinder_AoS");

    // SoA
    std::vector<mask<cylinder2D<false, soa::cylinder_intersector>,
                     std::uint_least16_t, algebra_v<scalar>>>
        cyl2_v;
    cyl2_v.emplace_back(0u, 0.5f * image.height() * rand, 0.5f * image.width(),
                        0.7f * image.width());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(cyl2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Cylinder SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "cylinder_SoA");

    // render a portal cylinder mask

    // AoS
    const mask<cylinder2D<false, cylinder_portal_intersector>,
               std::uint_least16_t, algebra_s<scalar>>
        pt_cyl2_s{0u, 0.1f * image.height(), 0.5f * image.width(),
                  0.7f * image.width()};
    std::vector<mask<cylinder2D<false, cylinder_portal_intersector>>>
        pt_cyl2_vec(soa_size, pt_cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_s>(image, std::move(pt_cyl2_vec),
                                           trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "\nPortal Cylinder AoS: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms" << std::endl;

    ppm.write(image, "portal_cylinder_AoS");

    // SoA
    std::vector<mask<cylinder2D<false, soa::cylinder_portal_intersector>,
                     std::uint_least16_t, algebra_v<scalar>>>
        pt_cyl2_v;
    pt_cyl2_v.emplace_back(0u, 0.5f * image.height() * rand,
                           0.5f * image.width(), 0.7f * image.width());

    start = std::chrono::high_resolution_clock::now();
    render_single_shape<scalar, algebra_v>(image, std::move(pt_cyl2_v), trfs_v,
                                           mat);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Portal Cylinder SoA: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                      start)
                         .count() /
                     1000'000.
              << " ms\n"
              << std::endl;

    ppm.write(image, "portal_cylinder_SoA");

    return EXIT_SUCCESS;
}
