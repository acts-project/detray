/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, read_detector) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;

    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string surface_file = data_directory + std::string("odd.csv");
    std::string layer_volume_file =
        data_directory + std::string("odd-layer-volumes.csv");
    std::string surface_grid_file =
        data_directory + std::string("odd-surface-grids.csv");
    std::string surface_grid_entries_file = "";

    auto d = detector_from_csv<>("odd", surface_file, layer_volume_file,
                                 surface_grid_file, surface_grid_entries_file,
                                 host_mr);

    std::cout << d.to_string() << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
