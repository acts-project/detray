/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"

#include <iostream>
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(__plugin, read_detector)
{
    using namespace detray;

    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string surface_file = data_directory + std::string("tml.csv");
    std::string surface_grid_file = data_directory + std::string("tml-surface-grids.csv");
    std::string layer_volume_file = data_directory + std::string("tml-layer-volumes.csv");

    auto d = detector_from_csv<static_transform_store>("tml", surface_file, surface_grid_file, layer_volume_file);

    std::cout << d.to_string() << std::endl;
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
