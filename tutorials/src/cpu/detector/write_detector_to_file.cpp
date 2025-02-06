// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

// Project include(s)
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/test/utils/detectors/build_toy_detector.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

/// Write a dector using the json IO
int main() {

    // First, create an example detector in in host memory to be written to disk
    vecmem::host_memory_resource host_mr;
    const auto [det, names] =
        detray::build_toy_detector<detray::tutorial::algebra_t>(host_mr);

    // Configuration for the writer:
    //     - use json format
    //     - replace the files if called multiple times
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);

    // Takes the detector 'det', a volume name map (only entry here the
    // detector name) and the writer config
    detray::io::write_detector(det, names, writer_cfg);
}
