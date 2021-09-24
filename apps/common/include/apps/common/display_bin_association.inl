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
#include "grids/associator.hpp"
#include "io/csv_io.hpp"
#include "style/styles.hpp"
#include "utils/enumerate.hpp"
#include "utils/generators.hpp"
#include "view/draw.hpp"
#include "view/views.hpp"

int main(int argc, char **argv) {
    vecmem::host_memory_resource host_mr;

    using point2 = __plugin::point2;
    using namespace detray;
    using namespace matplot;

    if (argc > 1) {
        std::string first_arg = argv[1];
        if (first_arg == "-h" or first_arg == "--help") {
            std::cout
                << "[detray] Usage: 'display_bin_association detector_name "
                   "<surface_file> <grid_file> <volume_file>'"
                << std::endl;
            return 1;
        } else if (argc > 7) {

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
            std::cout << "[detray] Detector read successfully." << std::endl;

            global_xy_view xy_view;
            global_z_phi_view zphi_view;

            // The layer, the bin
            int lvol = atoi(argv[6]);
            int bin_0 = atoi(argv[7]);
            int bin_1 = atoi(argv[8]);

            // view field additions:
            std::array<scalar, 4> grid_adds = {0., 0.2, 20., 0.};
            std::array<scalar, 4> cell_adds = {0., 0., 0., 0.};
            std::array<scalar, 4> assoc_adds = {5., 0.05, 5., 0.05};

            style surface_style;
            surface_style.fill_color = {0.9, 0.1, 0.6, 0.6};
            surface_style.line_width = 1;

            style selected_surface_style;
            selected_surface_style.fill_color = {0.2, 0.1, 0.6, 0.6};
            selected_surface_style.line_width = 1;

            style grid_style;
            grid_style.fill_color = {0., 0., 0., 0.};
            grid_style.line_width = 1;
            grid_style.line_style = ":";

            style cell_style;
            cell_style.fill_color = {0., 1.0, 0.0, 0.0};
            cell_style.line_width = 1;

            style assoc_style;
            assoc_style.fill_color = {0., 1.0, 0.0, 0.0};
            assoc_style.line_width = 1;
            assoc_style.line_style = "--";

            // Grid drawing sections styles and addons
            // - styles and addons are on faint grid, cell, assoc
            std::vector<style> gstyles = {grid_style, cell_style, assoc_style};
            std::vector<std::array<scalar, 4>> gadds = {grid_adds, cell_adds,
                                                        assoc_adds};

            decltype(d)::transform_store::context s_context;

            // Surface finders, volume, bounds
            auto surfaces_finders = d.surfaces_finders();
            const auto &lvolume = d.indexed_volume(lvol);
            const auto &bounds = lvolume.bounds();
            bool is_cylinder = std::abs(bounds[1] - bounds[0]) <
                               std::abs(bounds[3] - bounds[2]);
            const auto &sf_range = lvolume.template range<>();

            dindex finder_entry = lvolume.surfaces_finder_entry();
            const auto &surfaces = d.surfaces();
            const auto &surface_transforms = d.transforms(sf_range, s_context);
            const auto &surface_masks = d.masks();

            if (not is_cylinder) {

                const auto &disk_grid = surfaces_finders[finder_entry];
                auto r_borders = disk_grid.axis_p0().borders(bin_0);
                auto phi_borders = disk_grid.axis_p1().borders(bin_1);

                scalar r_min = r_borders[0];
                scalar r_max = r_borders[1];
                scalar phi_min = phi_borders[0];
                scalar phi_max = phi_borders[1];

                const auto &bin_entry = disk_grid.bin(bin_0, bin_1);
                std::cout << "Disk grid - bin (" << bin_0 << ", " << bin_1
                          << ") = ";
                for (auto b_e : bin_entry) {
                    std::cout << b_e << ", ";
                }
                std::cout << std::endl;

                // Loop over the surfaces within a volume
                for (dindex sfi = sf_range[0]; sfi < sf_range[1]; sfi++) {
                    const auto &s = surfaces[sfi];
                    dvector<point3> vertices = {};
                    const auto &mask_link = s.mask();
                    const auto &transform_link = s.transform();

                    // Unroll the mask container and generate vertices
                    const auto &transform = surface_transforms[transform_link];

                    const auto &mask_context = std::get<0>(mask_link);
                    const auto &mask_range = std::get<1>(mask_link);

                    auto vertices_per_masks = unroll_masks_for_vertices(
                        surface_masks, mask_range, mask_context,
                        std::make_integer_sequence<
                            dindex, std::tuple_size_v<decltype(
                                        d)::surface_mask_container>>{});

                    for (auto &vertices : vertices_per_masks) {
                        if (not vertices.empty()) {
                            style draw_style =
                                (std::find(
                                     bin_entry.begin(), bin_entry.end(),
                                     static_cast<dindex>(sfi - sf_range[0])) !=
                                 bin_entry.end())
                                    ? selected_surface_style
                                    : surface_style;

                            draw_vertices(vertices, transform, draw_style,
                                          xy_view);
                        }
                    }
                }

                // Grid drawing section
                for (auto [i, st] : enumerate(gstyles)) {

                    // Arc parameters
                    scalar arc_r_min = r_min - gadds[i][0];
                    scalar arc_r_max = r_max + gadds[i][0];
                    scalar arc_phi_min = phi_min - gadds[i][1];
                    scalar arc_phi_max = phi_max + gadds[i][1];
                    draw_arc(0., 0., arc_r_min, arc_phi_min, arc_phi_max, st);
                    draw_arc(0., 0., arc_r_max, arc_phi_min, arc_phi_max, st);

                    // Radial parameters
                    scalar rad_r_min = r_min - gadds[i][2];
                    scalar rad_r_max = r_max + gadds[i][2];
                    scalar rad_phi_min = phi_min - gadds[i][3];
                    scalar rad_phi_max = phi_max + gadds[i][3];
                    scalar cos_phi_min = std::cos(rad_phi_min);
                    scalar sin_phi_min = std::sin(rad_phi_min);
                    scalar cos_phi_max = std::cos(rad_phi_max);
                    scalar sin_phi_max = std::sin(rad_phi_max);

                    draw_line(rad_r_min * cos_phi_max, rad_r_min * sin_phi_max,
                              rad_r_max * cos_phi_max, rad_r_max * sin_phi_max,
                              st);
                    draw_line(rad_r_min * cos_phi_min, rad_r_min * sin_phi_min,
                              rad_r_max * cos_phi_min, rad_r_max * sin_phi_min,
                              st);
                }
            } else {
                const auto &cylinder_grid = surfaces_finders[finder_entry + 2];
                auto z_borders = cylinder_grid.axis_p0().borders(bin_0);
                auto phi_borders = cylinder_grid.axis_p1().borders(bin_1);

                scalar z_min = z_borders[0];
                scalar z_max = z_borders[1];
                scalar phi_min = phi_borders[0];
                scalar phi_max = phi_borders[1];

                const auto &bin_entry = cylinder_grid.bin(bin_0, bin_1);
                std::cout << "Cylinder grid - bin (" << bin_0 << ", " << bin_1
                          << ") = ";
                for (auto b_e : bin_entry) {
                    std::cout << b_e << ", ";
                }
                std::cout << std::endl;

                // Loop over the surfaces within a volume
                for (dindex sfi = sf_range[0]; sfi < sf_range[1]; sfi++) {
                    const auto &s = surfaces[sfi];
                    dvector<point3> vertices = {};
                    const auto &mask_link = s.mask();
                    const auto &transform_link = s.transform();

                    // Unroll the mask container and generate vertices
                    const auto &transform = surface_transforms[transform_link];

                    const auto &mask_context = std::get<0>(mask_link);
                    const auto &mask_range = std::get<1>(mask_link);

                    auto vertices_per_masks = unroll_masks_for_vertices(
                        surface_masks, mask_range, mask_context,
                        std::make_integer_sequence<
                            dindex, std::tuple_size_v<decltype(
                                        d)::surface_mask_container>>{});

                    for (auto &vertices : vertices_per_masks) {

                        if (not vertices.empty()) {

                            // Check the association (with potential splits)
                            style draw_style =
                                (std::find(
                                     bin_entry.begin(), bin_entry.end(),
                                     static_cast<dindex>(sfi - sf_range[0])) !=
                                 bin_entry.end())
                                    ? selected_surface_style
                                    : surface_style;

                            draw_vertices(vertices, transform, draw_style,
                                          zphi_view, true);
                        }
                    }
                }

                // Grid drawing section
                for (auto [i, st] : enumerate(gstyles)) {
                    if (i > 0) {
                        point2 p0 = {z_min - gadds[i][0],
                                     phi_min - gadds[i][1]};
                        point2 p1 = {z_min - gadds[i][0],
                                     phi_max + gadds[i][1]};
                        point2 p2 = {z_max + gadds[i][0],
                                     phi_max + gadds[i][1]};
                        point2 p3 = {z_max + gadds[i][0],
                                     phi_min - gadds[i][1]};

                        draw_polygon({p0, p1, p2, p3}, st);
                    }
                }
            }

            std::string vol_lay_name = "bin_assoc_";
            vol_lay_name += std::to_string(lvol);
            vol_lay_name += "_bin";
            vol_lay_name += std::to_string(bin_0);
            vol_lay_name += "_";
            vol_lay_name += std::to_string(bin_1);
            vol_lay_name += ".png";
            if (is_cylinder) {
                ax->xlabel("z [mm]");
                ax->ylabel("phi [rad]");
            } else {
                ax->xlabel("x [mm]");
                ax->ylabel("y [mm]");
                matplot::axis(equal);
            }

            ax->parent()->quiet_mode(false);
            show();

            save(vol_lay_name, true);

            return 1;
        }
    }

    std::cout << "[detray] Not enough arguments given, run with -h for help. "
              << std::endl;
    return 0;
}
