/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/proto_detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"

#include <iostream>
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(__plugin, read_detector)
{
    using namespace detray;

    std::string surface_file = "/Users/salzburg/Documents/work/dev/detray/data/tml.csv";
    std::string surface_grid_file = "/Users/salzburg/Documents/work/dev/detray/data/tml-surface-grids.csv";
    std::string layer_volume_file = "/Users/salzburg/Documents/work/dev/detray/data/tml-layer-volumes.csv";

    auto d = detector_from_csv<static_transform_store>("tml", surface_file, surface_grid_file, layer_volume_file);

}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

