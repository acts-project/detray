/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <iostream>
#include <map>
#include <string>

#include "core/detector.hpp"
#include "io/csv_io.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray_tests {

/** Keeps the relevant csv file names */
struct detector_input_files {
    std::string det_name, surface, layer_volume, surface_grid,
        surface_grid_entries;
};

/// open data detector
detector_input_files odd_files = {"odd", "odd.csv", "odd-layer-volumes.csv",
                                  "odd-surface-grids.csv", ""};

/// track ml detector
detector_input_files tml_files = {"tml", "tml.csv", "tml-layer-volumes.csv",
                                  "tml-surface-grids.csv", ""};

/** Read a detector from csv files */
auto read_from_csv(detector_input_files &files) {
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string surfaces = data_directory + files.surface;
    std::string volumes = data_directory + files.layer_volume;
    std::string grids = data_directory + files.surface_grid;
    std::string grid_entries = files.surface_grid_entries;
    std::map<detray::dindex, std::string> name_map{};

    auto d = detray::detector_from_csv<>(files.det_name, surfaces, volumes,
                                         grids, grid_entries, name_map);

    return std::make_pair<decltype(d), decltype(name_map)>(std::move(d),
                                                           std::move(name_map));
    ;
};

}  // namespace detray_tests