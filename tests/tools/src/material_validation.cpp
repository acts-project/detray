/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/track_generator_options.hpp"
#include "detray/test/detail/register_checks.hpp"
#include "detray/validation/detector_material_scan.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;
using namespace detray;

int main(int argc, char **argv) {

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<>;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    detray::material_scan<detector_t>::config mat_scan_cfg{};
    mat_scan_cfg.track_generator().uniform_eta(true);

    std::string description{"\ndetray material validation options"};
    detray::options::parse_options(description, argc, argv, reader_cfg,
                                   mat_scan_cfg.track_generator());

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    // Print the detector's material as recorded by a ray scan
    detray::detail::register_checks<detray::material_scan>(det, names,
                                                           mat_scan_cfg);

    // Run the checks
    return RUN_ALL_TESTS();
}
