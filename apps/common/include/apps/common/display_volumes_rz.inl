/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/proto_detector.hpp"
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
            auto d = detector_from_csv<static_transform_store>(name, surfaces, grids, volumes);
            std::cout << "[detray] Detector read successfully." << std::endl;
            std::cout << "         Volumes : " << d.volumes().size() << std::endl;
            // Parse the volumes for r/z max dimensions - pre-loop
            const scalar scalar_limit = std::numeric_limits<scalar>::max();
            darray<scalar, 4> rz_det = {scalar_limit, 0., scalar_limit, -scalar_limit};
            dmap<scalar, std::vector<dindex> > r_min_attachments;
            dmap<scalar, std::vector<dindex> > r_max_attachments;
            dmap<scalar, std::vector<dindex> > z_min_attachments;
            dmap<scalar, std::vector<dindex> > z_max_attachments;

            auto attach = [](dmap<scalar, std::vector<dindex> > &attachments, scalar value, dindex volume_index) -> void {
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

            auto ax = matplot::subplot({0.1, 0.1, 0.8, 0.8});
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

            // Helper method to sort and rmove duplicates
            auto sort_and_remove_duplicates = [](std::vector<dindex> &att) -> void {
                std::sort(att.begin(), att.end());
                att.erase(std::unique(att.begin(), att.end()), att.end());
            };

            std::cout << "[detray] Detector grid bins r_min " << r_min_attachments.size() << std::endl;
            std::cout << "                            r_max " << r_max_attachments.size() << std::endl;
            std::cout << "                            z_min " << z_min_attachments.size() << std::endl;
            std::cout << "                            z_max " << z_max_attachments.size() << std::endl;

            if (argc > 5)
            {
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

            show();

            return 1;
        }
    }

    std::cout << "[detray] Not enough arguments given, run with -h for help. " << std::endl;
    return 0;
}