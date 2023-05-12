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
#include "detray/masks/masks.hpp"

// Example linear algebra plugin: std::array
#include "detray/examples/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

/// Show how to build the test detectors.
///
/// Toy detector: Pixel section of the ACTS Generic detector (TrackML).
///
/// Telescope detector: Array of surfaces of a given type in a bounding portal
///                     box.
int main() {

    // Memory resource to allocate the detector data stores
    vecmem::host_memory_resource host_mr;

    //
    // Toy detector
    //

    // Toy detector type
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    // Number of barrel layers (0 - 4)
    constexpr std::size_t n_brl_layers{4};
    // Number of endcap layers on either side (0 - 7)
    // Note: The detector must be configured with 4 barrel layer to be able to
    // add encap layers
    constexpr std::size_t n_edc_layers{3};
    // Fill the detector
    const toy_detector_t toy_det =
        detray::create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Print the volume graph of the toy detector
    std::cout << detray::volume_graph{toy_det}.to_string() << std::endl;

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

    // Mask of a rectangle surfaces (20x20 mm)
    detray::mask<detray::rectangle2D<>> rectangle{
        0u, 20.f * detray::unit<detray::scalar>::mm,
        20.f * detray::unit<detray::scalar>::mm};

    // Build from given module positions
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

    // Case 1: Defaults: Straight telescope in z-direction,
    //         10 rectangle surfaces, 500mm in length, modules evenly spaced,
    //         no B-field
    const rectgl_telescope_t tel_det1 =
        detray::create_telescope_detector(host_mr, rectangle);
    // std::cout << detray::volume_graph{tel_det1}.to_string() << std::endl;

    // Case 2: Straight telescope in z-direction, 100 trapezoid surfaces, 2000mm
    //         in length, modules evenly spaced, no b-field
    const rectgl_telescope_t tel_det2 = detray::create_telescope_detector(
        host_mr, rectangle, 100, 2000.f * detray::unit<detray::scalar>::mm);
    // std::cout << detray::volume_graph{tel_det2}.to_string() << std::endl;
}
