/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/examples/types.hpp"
#include "detray/geometry/volume_graph.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

/// Show how to build the test detectors.
///
/// Toy Detector: Pixel section of the ACTS Generic detector (TrackML).
/// Telescope Detector: Array of surfaces of a given type in a bounding portal
///                     box.
int main() {

    //
    // Toy detector configuration
    //

    // Number of barrel layers (0 - 4)
    constexpr std::size_t n_brl_layers{4};
    // Number of endcap layers on either side (0 - 7)
    // Note: The detector must be configured with 4 barrel layer to be able to
    // add encap layers
    constexpr std::size_t n_edc_layers{3};
    // Memory resource to allocate the detector data stores
    vecmem::host_memory_resource host_mr;
    // Fill the detector
    auto toy_det =
        detray::create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Print the volume graph of the toy detector
    detray::volume_graph graph(toy_det);
    std::cout << graph.to_string() << std::endl;

    //
    // Telescope detector configuration
    //

    // Case 1:
}
