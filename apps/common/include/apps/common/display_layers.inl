/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <matplot/matplot.h>

#include <climits>
#include <fstream>
#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "style/styles.hpp"
#include "utils/enumerate.hpp"
#include "utils/generators.hpp"
#include "view/draw.hpp"
#include "view/views.hpp"

int main(int argc, char **argv) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace matplot;

    if (argc > 1) {
        std::string first_arg = argv[1];
        if (first_arg == "-h" or first_arg == "--help") {
            std::cout << "[detray] Usage: 'display_layers detector_name "
                         "<surface_file> <grid_file> <volume_file>'"
                      << std::endl;
            return 1;
        } else if (argc > 5) {

            int layer = argc > 6 ? atoi(argv[6]) : -1;
            bool mouse_over_mode = argc > 7 ? true : false;

            // Create a sub plot
            auto ax = matplot::subplot({0.1, 0.1, 0.65, 0.8});
            ax->parent()->quiet_mode(true);

            std::string name = first_arg;
            std::string surfaces_file = argv[2];
            std::string volumes_file = argv[3];
            std::string grids_file = argv[4];
            std::string grid_entries_file = argv[5];

            auto d =
                detector_from_csv<>(name, surfaces_file, volumes_file,
                                    grids_file, grid_entries_file, host_mr);

            std::cout << d.to_string() << std::endl;

            global_xy_view xy_view;
            global_z_phi_view zphi_view;

            style surface_style;
            surface_style.fill_color = {0.75, 0.1, 0.6, 0.6};
            surface_style.line_width = 1;

            style grid_style;
            grid_style.fill_color = {0., 1., 0., 0.};
            grid_style.line_width = 1;

            decltype(d)::transform_store::context s_context;

            auto surfaces_finders = d.surfaces_finders();

            // Loop over the volumes
            for (const auto [iv, v] : enumerate(d.volumes())) {

                if (layer > 0 and layer != iv) {
                    continue;
                }

                // Switch 'quiet' mode on during drawing
                const auto &bounds = v.bounds();
                bool is_cylinder = std::abs(bounds[1] - bounds[0]) <
                                   std::abs(bounds[3] - bounds[2]);
                const auto &sf_range = v.template range<>();

                const auto &surfaces = d.surfaces<>();
                const auto &volume_transforms =
                    d.transforms(sf_range, s_context);
                const auto &volume_masks = d.masks();

                if (v.empty()) {
                    continue;
                }

                // Loop over the surfaces within a volume
                for (dindex sfi = sf_range[0]; sfi < sf_range[1]; sfi++) {
                    const auto &s = surfaces[sfi];
                    dvector<point3> vertices = {};
                    const auto &mask_link = s.mask();
                    const auto &transform_link = s.transform();

                    // Unroll the mask container and generate vertices
                    const auto &transform = volume_transforms[transform_link];

                    const auto &mask_context = std::get<0>(mask_link);
                    const auto &mask_range = std::get<1>(mask_link);

                    auto vertices_per_masks = unroll_masks_for_vertices(
                        volume_masks, mask_range, mask_context,
                        std::make_integer_sequence<
                            dindex, std::tuple_size_v<decltype(
                                        d)::surface_mask_container>>{});

                    for (auto &vertices : vertices_per_masks) {
                        if (not vertices.empty()) {
                            if (not is_cylinder) {
                                draw_vertices(vertices, transform,
                                              surface_style, xy_view);
                            } else {
                                draw_vertices(vertices, transform,
                                              surface_style, zphi_view, true);
                            }
                        }
                    }
                }

                // Draw the surface finder grid
                dindex surfaces_finder_entry = v.surfaces_finder_entry();

                if (surfaces_finder_entry != dindex_invalid) {

                    if (not is_cylinder) {
                        const auto &r_phi_grid =
                            surfaces_finders[surfaces_finder_entry];
                        draw_r_phi_grid(r_phi_grid, grid_style);
                    } else {
                        const auto &z_phi_grid =
                            surfaces_finders[surfaces_finder_entry + 2];
                        draw_z_phi_grid(z_phi_grid, grid_style);
                    }
                }

                // Special functionality for single layers
                //
                if (layer > 0) {

                    // run option: mouse over parsing
                    if (mouse_over_mode) {
                        ax->run_command("set mouse; set print \"mouse.gpio\"");
                        ax->run_command("print MOUSE_X, MOUSE_Y, MOUSE_BUTTON");

                        std::string gpio_line;
                        std::ifstream gpio_file;
                        gpio_file.open("mouse.gpio", std::ios::in);
                        std::streamoff p = 0;

                        while (true) {
                            ax->run_command(
                                "print MOUSE_X, MOUSE_Y, MOUSE_BUTTON");

                            gpio_file.seekg(p);  //*1
                            while (getline(gpio_file, gpio_line)) {
                                std::cout << gpio_line << std::endl;

                                if (gpio_file.tellg() == -1)
                                    p = p + gpio_line.size();
                                else
                                    p = gpio_file.tellg();
                            }
                            gpio_file.clear();
                        }
                    }
                }
                std::string vol_lay_name = "lay_";
                vol_lay_name += std::to_string(iv);
                vol_lay_name += ".png";
                if (is_cylinder) {
                    ax->xlabel("z [mm]");
                    ax->ylabel("phi [rad]");
                } else {
                    ax->xlabel("x [mm]");
                    ax->ylabel("y [mm]");
                    matplot::axis(equal);
                }
                matplot::show();
                save(vol_lay_name, true);
            }

            return 1;
        }
    }

    std::cout << "[detray] Not enough arguments given, run with -h for help. "
              << std::endl;
    return 0;
}
