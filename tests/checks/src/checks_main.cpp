/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/checks/detector_consistency.hpp"
#include "detray/core/detector.hpp"
#include "detray/io/common/detector_reader.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <stdexcept>
#include <string>

template <typename check_t, typename detector_t>
void register_checks(const detector_t &det,
                     const typename detector_t::name_map &vol_names,
                     const std::string &test_name) {

    ::testing::RegisterTest("detray_check", test_name.c_str(), nullptr,
                            test_name.c_str(), __FILE__, __LINE__,
                            [=]() -> typename check_t::fixture_type * {
                                return new check_t(det, vol_names);
                            });
}

int main(int argc, char **argv) {

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Input data files
    detray::io::detector_reader_config reader_cfg{};
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            reader_cfg.add_file(argv[i]);
        }
    } else {
        throw std::runtime_error("Please specify an input file name!");
    }

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<>;

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    register_checks<detray::consistency_check<detector_t>>(
        det, names, "detector_consistency");

    // Run the checks
    return RUN_ALL_TESTS();
}
