/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Licenced under: Apache-2, see LICENSE file
 */

#include "core/surface.hpp"
#include "masks/rectangle2.hpp"
#include "utils/enumerate.hpp"

#include <fstream>
#include <map>
#include <string>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <vector>

/// @note plugin has to be defined with a preprocessor command

namespace detray
{

    template <typename surface_type>
    struct csv_layer
    {
        std::vector<surface_type> surfaces;
        std::vector<rectangle2<scalar>> rectangle_masks;
        //std::vector<trapezoid2> trapezoid_masks;
    };

    template <typename surface_type>
    struct csv_volume
    {
        std::vector<csv_layer<surface_type>> layers;
    };

    template <typename surface_type>
    struct csv_detector
    {
        std::vector<csv_volume<surface_type>> volumes;
    };

    template <typename transform_type>
    csv_detector<surface<transform_type, unsigned int, unsigned int>> read_csv(std::string filename)
    {

        using context = typename transform_type::context;
        using surface = surface<transform_type, unsigned int, unsigned int>;
        using rectangle2 = rectangle2<scalar>;

        struct surface_info
        {
            transform_type transform_info;
            std::array<scalar, 3> mask_info;
            unsigned int source_info;
        };

        std::map<unsigned int, std::map<unsigned int, std::map<unsigned int, surface_info>>> read_surface_info;
        std::ifstream csv_file(filename);

        if (!csv_file.is_open())
        {
            throw std::runtime_error("Could not open file");
        }

        std::string line;
        std::stringstream ss(line);
        // First line is header
        std::getline(csv_file, line);

        csv_detector<surface> detector;
        while (std::getline(csv_file, line))
        {
            std::vector<std::string> split_strings;
            std::stringstream ss(line);
            while (ss.good())
            {
                std::string sub_string;
                std::getline(ss, sub_string, ',');
                split_strings.push_back(sub_string);
            }
            unsigned int geo_id = static_cast<unsigned int>(std::atoi(split_strings[0].c_str()));
            int volume_id = std::atoi(split_strings[1].c_str());
            int layer_id = std::atoi(split_strings[2].c_str());
            int module_id = std::atoi(split_strings[3].c_str());

            float t0 = std::atof(split_strings[4].c_str());
            float t1 = std::atof(split_strings[5].c_str());
            float t2 = std::atof(split_strings[6].c_str());
            float x0 = std::atof(split_strings[7].c_str());
            float x1 = std::atof(split_strings[8].c_str());
            float x2 = std::atof(split_strings[9].c_str());
            float y0 = std::atof(split_strings[10].c_str());
            float y1 = std::atof(split_strings[11].c_str());
            float y2 = std::atof(split_strings[12].c_str());
            float z0 = std::atof(split_strings[13].c_str());
            float z1 = std::atof(split_strings[14].c_str());
            float z2 = std::atof(split_strings[15].c_str());

            std::array<scalar, 16> matrix_array = {x0, y0, z0, t0, x1, y1, z1, t1, x2, y2, z2, t2, 0., 0., 0., 1.};

            float mhx_miny = std::atof(split_strings[17].c_str());
            float mhx_maxy = std::atof(split_strings[18].c_str());
            float mhy = std::atof(split_strings[19].c_str());
            std::array<scalar, 3> mask = {mhx_miny, mhx_maxy, mhy};

            context default_context;
            read_surface_info[volume_id][layer_id].insert({module_id, {transform_type{matrix_array, default_context}, mask, geo_id}});
        }
        csv_file.close();

        // Now create a csv_detector representation
        for (auto [vkey, vvalue] : read_surface_info)
        {
            csv_volume<surface> volume;
            for (auto [lkey, lvalue] : vvalue)
            {
                std::vector<rectangle2> rectangle_masks;

                auto find_mask_index = [](std::vector<rectangle2> &masks, const std::array<scalar, 2> &mask_input) -> size_t {
                    for (auto [ i, rect ] : enumerate(masks)){
                        if (rect == mask_input){
                            return i;
                        }
                    }
                    auto i = masks.size();
                    rectangle2 r = { mask_input };
                    masks.push_back(r);
                    return i;
                };

                csv_layer<surface> layer;
                for (auto [mkey, mvalue] : lvalue)
                {
                    auto maks_index = find_mask_index(rectangle_masks, { mvalue.mask_info[0], mvalue.mask_info[2]} );
                    layer.surfaces.push_back(surface(std::move(mvalue.transform_info), std::move(maks_index), std::move(mvalue.source_info)));
                }
                layer.rectangle_masks = rectangle_masks;
                volume.layers.push_back(std::move(layer));
            }
            detector.volumes.push_back(volume);
        }
        return detector;
    }

} // namespace detray
