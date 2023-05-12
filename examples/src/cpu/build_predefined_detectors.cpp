/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracks/tracks.hpp"

// Example linear algebra plugin: std::array
#include "detray/examples/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <cstdlib>
#include <iostream>

/// Show how to build the test detectors.
///
/// Toy detector: Pixel section of the ACTS Generic detector (TrackML).
///
/// Telescope detector: Array of surfaces of a given type in a bounding portal
///                     box.
int main(int argc, char** argv) {

    // Memory resource to allocate the detector data stores
    vecmem::host_memory_resource host_mr;

    //
    // Toy detector
    //

    // Toy detector type
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    // Number of barrel layers (0 - 4)
    std::size_t n_brl_layers{4u};
    // Number of endcap layers on either side (0 - 7)
    // Note: The detector must be configured with 4 barrel layers to be able to
    // add any encap layers
    std::size_t n_edc_layers{1u};

    // Read toy detector config from commandline, if it was given
    if (argc == 3) {
        n_brl_layers = std::atoi(argv[1]);
        n_edc_layers = std::atoi(argv[2]);
    }

    // Fill the detector
    const toy_detector_t toy_det =
        detray::create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Print the volume graph of the toy detector
    std::cout << "\nToy detector:\n"
              << "-------------\n"
              << detray::volume_graph{toy_det}.to_string() << std::endl;

    //
    // Telescope detector
    //

    // The telescope detector is built according to a 'pilot trajectory'
    // (either a helix or a ray) along which the modules are placed in given
    // distances. the world portals are constructed from a bounding box around
    // the test surfaces (the envelope is configurable)

    // Telescope detector tyoes containing the given module shape (can be any)
    using rectgl_telescope_t =
        detray::detector<detray::telescope_metadata<detray::rectangle2D<>>>;
    using trapzd_telescope_t =
        detray::detector<detray::telescope_metadata<detray::trapezoid2D<>>>;

    // Mask with a rectangular shape (20x20 mm)
    detray::mask<detray::rectangle2D<>> rectangle{
        0u, 20.f * detray::unit<detray::scalar>::mm,
        20.f * detray::unit<detray::scalar>::mm};

    // Mask with a trapezoid shape
    constexpr detray::scalar hx_min_y{10.f * detray::unit<detray::scalar>::mm};
    constexpr detray::scalar hx_max_y{30.f * detray::unit<detray::scalar>::mm};
    constexpr detray::scalar hy{20.f * detray::unit<detray::scalar>::mm};
    constexpr detray::scalar divisor{10.f / (20.f * hy)};
    detray::mask<detray::trapezoid2D<>> trapezoid{0u, hx_min_y, hx_max_y, hy,
                                                  divisor};

    // Build from given module positions (places 11 surfaces)
    std::vector<detray::scalar> positions = {
        0.f * detray::unit<detray::scalar>::mm,
        50.f * detray::unit<detray::scalar>::mm,
        100.f * detray::unit<detray::scalar>::mm,
        150.f * detray::unit<detray::scalar>::mm,
        200.f * detray::unit<detray::scalar>::mm,
        250.f * detray::unit<detray::scalar>::mm,
        300.f * detray::unit<detray::scalar>::mm,
        350.f * detray::unit<detray::scalar>::mm,
        400.f * detray::unit<detray::scalar>::mm,
        450.f * detray::unit<detray::scalar>::mm,
        500.f * detray::unit<detray::scalar>::mm};

    //
    // Case 1: Defaults: Straight telescope in z-direction,
    //         10 rectangle surfaces, 500mm in length, modules evenly spaced,
    //         no B-field, silicon material (80mm)
    const rectgl_telescope_t tel_det1 =
        detray::create_telescope_detector(host_mr, rectangle);

    std::cout << "\nTelescope detector - case 1:\n"
              << "----------------------------\n"
              << detray::volume_graph{tel_det1}.to_string() << std::endl;

    //
    // Case 2: Straight telescope in z-direction, 15 trapezoid surfaces, 2000mm
    //         in length, modules evenly spaced, no B-field,
    //         silicon material (80mm)
    const trapzd_telescope_t tel_det2 = detray::create_telescope_detector(
        host_mr, trapezoid, 15, 2000.f * detray::unit<detray::scalar>::mm);

    std::cout << "\nTelescope detector - case 2:\n"
              << "----------------------------\n"
              << detray::volume_graph{tel_det2}.to_string() << std::endl;

    //
    // Case 3: Straight telescope in x-direction, 11 rectangle surfaces, 2000mm
    //         in length, modules evenly spaced, no B-field,
    //         silicon material (80mm)

    // Pilot trajectory in x-direction
    detray::detail::ray<detray::example::transform3> x_track{
        {0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};

    const rectgl_telescope_t tel_det3 = create_telescope_detector(
        host_mr, rectangle, positions, detray::silicon_tml<detray::scalar>(),
        80.f * detray::unit<detray::scalar>::um, x_track);

    std::cout << "\nTelescope detector - case 3:\n"
              << "----------------------------\n"
              << detray::volume_graph{tel_det3}.to_string() << std::endl;

    //
    // Case 4: Bent telescope along helical track, 11 trapezoid surfaces,
    //         modules spaced according to given positions,
    //         constant B-field, silicon material (80mm)

    // Pilot track in x-direction
    detray::free_track_parameters<detray::example::transform3> y_track{
        {0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
    // Helix in a constant B-field in z-direction
    detray::example::vector3 B_z{0.f, 0.f,
                                 1.f * detray::unit<detray::scalar>::T};
    detray::detail::helix<detray::example::transform3> helix(y_track, &B_z);

    // Prepare constant B-field
    using b_field_t = typename trapzd_telescope_t::bfield_type;
    b_field_t b_field_z{
        b_field_t::backend_t::configuration_t{B_z[0], B_z[1], B_z[2]}};

    const trapzd_telescope_t tel_det4 = detray::create_telescope_detector(
        host_mr, std::move(b_field_z), trapezoid, positions,
        detray::silicon_tml<detray::scalar>(),
        80.f * detray::unit<detray::scalar>::um, helix);

    std::cout << "\nTelescope detector - case 4:\n"
              << "----------------------------\n"
              << detray::volume_graph{tel_det4}.to_string() << std::endl;
}
