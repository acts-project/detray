/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include <boost/program_options.hpp>

// System include(s)
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;
using namespace detray;

int main(int argc, char **argv) {

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<>;

    // Visualization style to be applied to the svgs
    auto style = detray::svgtools::styling::tableau_colorblind;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Options parsing
    po::options_description desc("\ndetray detector validation options");

    std::vector<dindex> volumes, surfaces;
    desc.add_options()("help", "produce help message")(
        "outdir", po::value<std::string>(), "Output directory for plots")(
        "geometry_file", po::value<std::string>(), "Geometry input file")(
        "volumes", po::value<std::vector<dindex>>(&volumes)->multitoken(),
        "List of volumes that should be displayed")(
        "surfaces", po::value<std::vector<dindex>>(&surfaces)->multitoken(),
        "List of surfaces that should be displayed")(
        "hide_portals", "Hide portal surfaces")("hide_passives",
                                                "Hide passive surfaces")(
        "show_info", "Show info boxes")("write_volume_graph",
                                        "Writes the volume graph to file");

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc,
                                 po::command_line_style::unix_style ^
                                     po::command_line_style::allow_short),
              vm);
    po::notify(vm);

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};

    // Help message
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    // General options
    std::string outdir{vm.count("outdir") ? vm["outdir"].as<std::string>()
                                          : "./plots/"};
    auto path = std::filesystem::path(outdir);
    if (not std::filesystem::exists(path)) {
        std::error_code err;
        if (!std::filesystem::create_directories(path, err)) {
            throw std::runtime_error(err.message());
        }
    }

    // Input files
    if (vm.count("geometry_file")) {
        reader_cfg.add_file(vm["geometry_file"].as<std::string>());
    } else {
        std::stringstream err_stream{};
        err_stream << "Please specify a geometry input file!\n\n" << desc;

        throw std::invalid_argument(err_stream.str());
    }

    // Read the detector geometry
    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);
    detector_t::geometry_context gctx{};

    // Creating the svg generator for the detector.
    detray::svgtools::illustrator il{det, names, style};
    il.show_info(vm.count("show_info"));
    il.hide_portals(vm.count("hide_portals"));
    il.hide_passives(vm.count("hide_passives"));

    actsvg::style::stroke stroke_black = actsvg::style::stroke();

    // x-y axis.
    auto xy_axis = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                          stroke_black, "x", "y");
    // z-r axis.
    auto zr_axis = actsvg::draw::x_y_axes("axes", {-2000, 2000}, {-5, 250},
                                          stroke_black, "z", "r");

    // Creating the views.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // Display the volumes
    if (not volumes.empty()) {
        std::string name_xy = names.at(0) + "_volume_display_xy";
        const auto vol_xy_svg = il.draw_volumes(name_xy, volumes, xy, gctx);
        detray::svgtools::write_svg(path / name_xy, {xy_axis, vol_xy_svg});

        std::string name_zr = names.at(0) + "_volume_display_zr";
        const auto vol_zr_svg = il.draw_volumes(name_zr, volumes, zr, gctx);
        detray::svgtools::write_svg(path / name_zr, {zr_axis, vol_zr_svg});
    }

    // Display the surfaces
    if (not surfaces.empty()) {
        std::string name_xy = names.at(0) + "_surface_display_xy";
        const auto sf_xy_svg = il.draw_surfaces(name_xy, surfaces, xy, gctx);
        detray::svgtools::write_svg(name_xy, {xy_axis, sf_xy_svg});

        std::string name_zr = names.at(0) + "_surface_display_zr";
        const auto sf_zr_svg = il.draw_surfaces(name_zr, surfaces, zr, gctx);
        detray::svgtools::write_svg(path / name_zr, {zr_axis, sf_zr_svg});
    }

    // If nothing was specified, display the whole detector
    if (volumes.empty() and surfaces.empty()) {
        std::string name_xy = names.at(0) + "_display_xy";
        const auto det_xy_svg = il.draw_detector(name_xy, xy, gctx);
        detray::svgtools::write_svg(path / name_xy, {xy_axis, det_xy_svg});

        std::string name_zr = names.at(0) + "_display_zr";
        const auto det_zr_svg = il.draw_detector(name_zr, zr, gctx);
        detray::svgtools::write_svg(path / name_zr, {zr_axis, det_zr_svg});
    }

    // Display the detector volume graph
    if (vm.count("write_volume_graph")) {
        detray::volume_graph graph(det);

        detray::io::detail::file_handle stream{
            path / (names.at(0) + "_volume_graph.dot"),
            std::ios::out | std::ios::trunc};
        *stream << graph.to_dot_string();
    }
}
