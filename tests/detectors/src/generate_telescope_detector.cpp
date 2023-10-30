/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/io/common/detector_writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Boost
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace detray;

namespace {

/// Generate and write a telescope detector, given the commandline variables
/// and a configuration for the detector writer @param writer_cfg
template <typename mask_shape_t, typename value_t>
void write_telecope(const po::variables_map &vm,
                    io::detector_writer_config &writer_cfg,
                    std::vector<value_t> &mask_params) {

    using detector_t = detector<telescope_metadata<mask_shape_t>>;
    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;

    detray::tel_det_config<mask_shape_t> tel_cfg{mask_params};

    tel_cfg.n_surfaces(vm["modules"].as<unsigned int>());
    tel_cfg.length(vm["length"].as<scalar_t>());
    tel_cfg.mat_thickness(vm["thickness"].as<scalar_t>());

    std::string direction{vm["direction"].as<std::string>()};
    detray::detail::ray<transform3_t> r{};
    if (direction == "x") {
        r.set_dir({1.f, 0.f, 0.f});
    } else if (direction == "y") {
        r.set_dir({0.f, 1.f, 0.f});
    } else if (direction == "z") {
        r.set_dir({0.f, 0.f, 1.f});
    }
    tel_cfg.pilot_track(r);

    // Build the detector
    vecmem::host_memory_resource host_mr;
    auto [tel_det, tel_names] = create_telescope_detector(host_mr, tel_cfg);

    // Write to file
    detray::io::write_detector(tel_det, tel_names, writer_cfg);
}

}  // anonymous namespace

int main(int argc, char **argv) {

    using scalar_t = detray::scalar;

    // Options parsing
    po::options_description desc("\nTelescope detector generation options");

    std::vector<scalar_t> mask_params{
        20.f * unit<scalar_t>::mm,
        20.f * unit<scalar_t>::mm};  // < default values for rectangles
    desc.add_options()("help", "produce help message")(
        "outdir", po::value<std::string>(), "Output directory for files")(
        "write_volume_graph", "writes the volume graph to file")(
        "compactify_json", "not implemented")(
        "write_material", "toggle material output")("write_grids",
                                                    "toggle grid output")(
        "modules", po::value<unsigned int>()->default_value(10u),
        "number of modules in telescope [1-20]")(
        "type", po::value<std::string>()->default_value("rectangle"),
        "type of the telescope modules [rectangle, trapezoid, annulus, ring, "
        "cylinder]")(
        "params", po::value<std::vector<scalar_t>>(&mask_params)->multitoken(),
        "Mask values for the shape given in 'type'")(
        "length",
        po::value<scalar_t>()->default_value(500.f * unit<scalar_t>::mm),
        "length of the telescope [mm]")(
        "thickness",
        po::value<scalar_t>()->default_value(1.f * unit<scalar_t>::mm),
        "thickness of the module silicon material")(
        "direction", po::value<std::string>()->default_value("z"),
        "direction of the telescope in global frame [x, y, z]");

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc,
                                 po::command_line_style::unix_style ^
                                     po::command_line_style::allow_short),
              vm);
    po::notify(vm);

    // Help message
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    // Configuration
    detray::io::detector_writer_config writer_cfg{};
    writer_cfg.format(detray::io::format::json).replace_files(false);

    // General options
    std::string outdir{vm.count("outdir") ? vm["outdir"].as<std::string>()
                                          : "./telescope_detector/"};
    writer_cfg.path(std::move(outdir));
    writer_cfg.compactify_json(vm.count("compactify_json"));
    writer_cfg.write_material(vm.count("write_material"));
    writer_cfg.write_grids(vm.count("write_grids"));

    // Build the geometry
    std::string type{vm["type"].as<std::string>()};
    if (type == "rectangle") {
        write_telecope<rectangle2D<>>(vm, writer_cfg, mask_params);
    } else if (type == "trapezoid") {
        write_telecope<trapezoid2D<>>(vm, writer_cfg, mask_params);
    } else if (type == "annulus") {
        write_telecope<annulus2D<>>(vm, writer_cfg, mask_params);
    } else if (type == "ring") {
        write_telecope<ring2D<>>(vm, writer_cfg, mask_params);
    } else if (type == "cylinder") {
        write_telecope<cylinder2D<>>(vm, writer_cfg, mask_params);
    }
}
