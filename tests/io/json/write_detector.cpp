/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "plugins/array_defs.hpp"
#include "tests/common/single_layer_detector.hpp"
#include "io/json_core.hpp"

#include <fstream>

int main(int argc, char **argv)
{
    auto d = createDetector();

    std::ofstream output_file;
    output_file.open(d.name()+std::string(".json"));
    output_file << json::write_detector(d).dump(2);

    output_file.close();

}
