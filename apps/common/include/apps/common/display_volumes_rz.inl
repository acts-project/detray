/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "utils/containers.hpp"

#include <iostream>
#include <climits>
#include <matplot/matplot.h>

int main(int argc, char **argv)
{
    using namespace detray;
    using namespace matplot;

    if (argc > 1)
    {
        std::string first_arg = argv[1];
        if (first_arg == "-h" or first_arg == "--help")
        {
            std::cout << "[detray] Usage: 'display_volumes_rz detector_name <surface_file> <grid_file> <volume_file>'" << std::endl;
            return 1;
        }
        else if (argc > 4)
        {
            std::string name = first_arg;
            std::string surfaces = argv[2];
            std::string grids = argv[3];
            std::string volumes = argv[4];
            auto d = detector_from_csv<static_transform_store<>>(name, surfaces, grids, volumes);
            std::cout << "[detray] Detector read successfully." << std::endl;
            std::cout << d.to_string() << std::endl;
            // std::cout << "         Volumes : " << d.volumes().size() << std::endl;
            // Parse the volumes for r/z max dimensions - pre-loop
            const scalar scalar_limit = std::numeric_limits<scalar>::max();
            darray<scalar, 4> rz_det = {scalar_limit, 0., scalar_limit, -scalar_limit};
            dmap<scalar, std::vector<dindex>> r_min_attachments;
            dmap<scalar, std::vector<dindex>> r_max_attachments;
            dmap<scalar, std::vector<dindex>> z_min_attachments;
            dmap<scalar, std::vector<dindex>> z_max_attachments;

            auto attach = [](dmap<scalar, std::vector<dindex>> &attachments, scalar value, dindex volume_index) -> void
            {
                if (attachments.find(value) == attachments.end())
                {
                    attachments[value] = {volume_index};
                }
                else
                {
                    attachments[value].push_back(volume_index);
                }
            };

            for (const auto &v : d.volumes())
            {
                const auto &v_bounds = v.bounds();
                rz_det[0] = std::min(rz_det[0], v_bounds[0]);
                rz_det[1] = std::max(rz_det[1], v_bounds[1]);
                rz_det[2] = std::min(rz_det[2], v_bounds[2]);
                rz_det[3] = std::max(rz_det[3], v_bounds[3]);
                attach(r_min_attachments, v_bounds[0], v.index());
                attach(r_max_attachments, v_bounds[1], v.index());
                attach(z_min_attachments, v_bounds[2], v.index());
                attach(z_max_attachments, v_bounds[3], v.index());
            }

            // Helper method to sort and rmove duplicates
            auto sort_and_remove_duplicates = [](std::vector<dindex> &att) -> void
            {
                std::sort(att.begin(), att.end());
                att.erase(std::unique(att.begin(), att.end()), att.end());
            };

            std::cout << "[detray] Detector grid bins r_min " << r_min_attachments.size() << std::endl;
            std::cout << "                            r_max " << r_max_attachments.size() << std::endl;
            std::cout << "                            z_min " << z_min_attachments.size() << std::endl;
            std::cout << "                            z_max " << z_max_attachments.size() << std::endl;

            auto ax = matplot::subplot({0.1, 0.1, 0.8, 0.8});
            ax->parent()->quiet_mode(true);

            if (argc > 5)
            {
                int base_draw_option = std::atoi(argv[5]);

                // Draw first the detailed views
                if (argc > 6)
                {
                    // drawing: base draw option
                    int detail_draw_option = std::atoi(argv[6]);

                    if (detail_draw_option == 1)
                    {
                        std::cout << "[detray] Drawing grid lines." << std::endl;

                        // Drawing the lines for the grid search
                        for (auto [key, value] : r_min_attachments)
                        {
                            sort_and_remove_duplicates(value);
                            auto l = line(rz_det[2], key, rz_det[3], key);
                            l->color("red");
                        }
                        for (auto [key, value] : r_max_attachments)
                        {
                            sort_and_remove_duplicates(value);
                            auto l = line(rz_det[2], key, rz_det[3], key);
                            l->color("red");
                        }
                        for (auto [key, value] : z_min_attachments)
                        {
                            sort_and_remove_duplicates(value);
                            auto l = line(key, rz_det[0], key, rz_det[1]);
                            l->color("red");
                        }
                        for (auto [key, value] : z_max_attachments)
                        {
                            sort_and_remove_duplicates(value);
                            auto l = line(key, rz_det[0], key, rz_det[1]);
                            l->color("red");
                        }
                    }

                    if (detail_draw_option == 2 and base_draw_option == 0)
                    {
                        std::cout << "[detray] Drawing volume heat map." << std::endl;

                        // Get the grid
                        auto volume_grid = d.volume_search_grid();

                        const auto &r_axis = volume_grid.axis_p0();
                        const auto &z_axis = volume_grid.axis_p1();
                        dvector<dvector<dindex>> heat_data;

                        heat_data.reserve(r_axis.bins());

                        for (dindex ir = 0; ir < r_axis.bins(); ++ir)
                        {
                            dvector<dindex> heat_values;
                            heat_values.reserve(z_axis.bins());
                            for (dindex iz = 0; iz < z_axis.bins(); ++iz)
                            {
                                dindex vidx = volume_grid.bin(ir, iz);
                                if (vidx % 2)
                                {
                                    vidx += 100;
                                }
                                heat_values.push_back(vidx);
                            }
                            heat_data.push_back(heat_values);
                        }
                        heatmap(heat_data);
                    }

                    // Get the grid
                    auto volume_grid = d.volume_search_grid();
                    const auto &r_axis = volume_grid.axis_p0();
                    const auto &z_axis = volume_grid.axis_p1();

                    // Volume grid lines
                    if (detail_draw_option == 3)
                    {
                        std::cout << "[detray] Drawing volume grid." << std::endl;

                        for (dindex ir = 0; ir < r_axis.bins(); ++ir)
                        {
                            darray<scalar, 2> r_borders = r_axis.borders(ir);
                            for (dindex iz = 0; iz < z_axis.bins(); ++iz)
                            {
                                darray<scalar, 2> z_borders = z_axis.borders(iz);
                                dindex volume = volume_grid.bin(ir, iz);

                                auto v_dis = rectangle(z_borders[0], r_borders[0], z_borders[1] - z_borders[0], r_borders[1] - r_borders[0], 0);
                                v_dis->color({0.9, 0.5, 0.5, 0.5});
                            }
                        }
                    }

                    // Volume grid colored
                    if (detail_draw_option == 4)
                    {
                        std::cout << "[detray] Drawing grid colored" << std::endl;
                        rectangle(rz_det[2], rz_det[0], rz_det[3] - rz_det[2], rz_det[1] - rz_det[0], 0);

                        std::map<dindex, darray<float, 3>> volume_colors;

                        /// Helper method for randomized color
                        auto random_color = []() -> darray<float, 3>
                        {
                            int r = std::rand() % 10;
                            int g = std::rand() % 10;
                            int b = std::rand() % 10;
                            return {r * 0.1f, g * 0.1f, b * 0.1f};
                        };

                        for (dindex ir = 0; ir < r_axis.bins(); ++ir)
                        {
                            for (dindex iz = 0; iz < z_axis.bins(); ++iz)
                            {
                                dindex volume_index = volume_grid.bin(ir, iz);
                                if (volume_colors.find(volume_index) == volume_colors.end())
                                {
                                    volume_colors[volume_index] = random_color();
                                }
                                auto r_borders = r_axis.borders(ir);
                                auto z_borders = z_axis.borders(iz);
                                auto rz = rectangle(z_borders[0], r_borders[0], z_borders[1] - z_borders[0], r_borders[1] - r_borders[0]);
                                auto color = volume_colors[volume_index];
                                rz->fill(true);
                                rz->color(color);
                            }
                        }
                    }

                    // Volume grid blocks
                    if (detail_draw_option == 5)
                    {

                        std::cout << "[detray] Drawing expanded volume grid around kernels" << std::endl;
                        rectangle(rz_det[2], rz_det[0], rz_det[3] - rz_det[2], rz_det[1] - rz_det[0], 0);

                        using bin2 = darray<dindex, 2>;
                        using bin_area = darray<bin2, 2>;

                        /** Helper method for expanding the area 
                         * 
                         * @param be The bin area to be expanded
                         * @param bin The bin in question  
                         */
                        auto expand = [&](bin_area &be, const bin2 &bin) -> void
                        {
                            auto &b = be[0];
                            auto &e = be[1];
                            b[0] = std::min(b[0], bin[0]);
                            b[1] = std::min(b[1], bin[1]);
                            e[0] = std::max(e[0], bin[0]);
                            e[1] = std::max(e[1], bin[1]);
                        };

                        bin2 max_bin = {r_axis.bins(), z_axis.bins()};
                        bin2 min_bin = {0, 0};
                        bin_area start_bin = {max_bin, min_bin};

                        dmap<dindex, bin_area> volume_extends;
                        for (dindex ir = 0; ir < r_axis.bins(); ++ir)
                        {
                            for (dindex iz = 0; iz < z_axis.bins(); ++iz)
                            {
                                dindex volume_index = volume_grid.bin(ir, iz);
                                auto vol_itr = volume_extends.find(volume_index);
                                if (vol_itr == volume_extends.end())
                                {
                                    volume_extends[volume_index] = start_bin;
                                }
                                expand(volume_extends[volume_index], {ir, iz});
                            }
                        }

                        std::cout << "[detray] Found " << volume_extends.size() << " individual volume extends." << std::endl;
                        for (const auto &[vidx, be] : volume_extends)
                        {
                            const auto &b = be[0];
                            const auto &e = be[1];
                            scalar vol_rmin = r_axis.borders(b[0])[0];
                            scalar vol_zmin = z_axis.borders(b[1])[0];
                            scalar vol_rmax = r_axis.borders(e[0])[1];
                            scalar vol_zmax = z_axis.borders(e[1])[1];
                            std::cout << "     Volume " << vidx << " with rmin/rmax = " << vol_rmin << ", " << vol_rmax << std::endl;
                            std::cout << "            " << vidx << " with zmin/zmax = " << vol_zmin << ", " << vol_zmax << std::endl;
                            std::cout << "            " << vidx << " from  = (" << b[0] << ", " << b[1] << ")" << std::endl;
                            std::cout << "            " << vidx << " to    = (" << e[0] << ", " << e[1] << ")" << std::endl;
                            auto vis_vol = rectangle(vol_zmin, vol_rmin, vol_zmax - vol_zmin, vol_rmax - vol_rmin, 0);
                        }
                    }
                }

                // Drawing: base draw option on top
                if (base_draw_option == 1)
                {

                    std::cout << "[detray] Drawing volumes." << std::endl;
                    rectangle(rz_det[2], rz_det[0], rz_det[3] - rz_det[2], rz_det[1] - rz_det[0], 0);

                    // Drawing loop for all volumes
                    for (const auto &v : d.volumes())
                    {
                        const auto &v_bounds = v.bounds();
                        auto v_dis = rectangle(v_bounds[2], v_bounds[0], v_bounds[3] - v_bounds[2], v_bounds[1] - v_bounds[0], 0);
                        v_dis->fill(true);
                        if (not v.empty())
                        {
                            v_dis->color({0.8, 0, 0, 1});
                        }
                        else
                        {
                            v_dis->color({0.9, 0.5, 0.5, 0.5});
                        }
                    }
                }
            }

            ax->parent()->quiet_mode(false);
            show();

            return 1;
        }
    }

    std::cout << "[detray] Not enough arguments given, run with -h for help. " << std::endl;
    return 0;
}