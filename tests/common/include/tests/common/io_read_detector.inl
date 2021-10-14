/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <string>

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, read_detector) {
    using namespace detray;

    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    /*std::string name = "odd";
    std::string surfaces = data_directory + "odd.csv";
    std::string volumes = data_directory + "odd-layer-volumes.csv";
    std::string grids = data_directory + "odd-surface-grids.csv";
    std::string grid_entries = "";*/

    std::string name = "tml";
    std::string surfaces = data_directory + "tml.csv";
    std::string volumes = data_directory + "tml-layer-volumes.csv";
    std::string grids = data_directory + "tml-surface-grids.csv";
    std::string grid_entries = "";
    std::map<dindex, std::string> name_map{};

    auto d = detector_from_csv<>(name, surfaces, volumes, grids, grid_entries,
                                 name_map);

    std::cout << d.to_string(name_map) << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
