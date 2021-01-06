/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "json/json_masks.hpp"

#include <algorithm>

#include <nlohmann/json.hpp>

using njson = nlohmann::json;

namespace detray
{

    namespace json
    {

        /** Adding an axis to the json output format
         * 
         * @tparam mask_type is the type of the axis object and @param a the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename axis_type>
        njson write_axis(const axis_type &a)
        {
            njson aj;
            aj["bins"] = a.axis_bins;
            aj["identifier"] = a.axis_identifier;
            aj["range"] = a.range();
            return aj;
        }

        /** Write the typed mask container to json using variadic templates
         * 
         * @tparam axis_type the type of the tuple entries
         * @param axes the axisd tuple
         * 
         * @return a valid json object with substructure
         **/
        template <typename... axis_type>
        njson write_all_axes(std::tuple<axis_type...> const &axes)
        {
            njson aj;
            std::apply(
                [&](axis_type const &... a) {
                    ((aj.push_back(write_axis(a)), ...));
                },
                axes);
            return aj;
        }

        /** Adding a grid to the json output format
         * 
         * @tparam grid_type is the type of the grid object and @param g the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename grid_type>
        njson write_grid(const grid_type &g)
        {
            njson gj;

            auto axes = g.axes();
            gj["axes"] = write_all_axes(axes);
            return gj;
        }

    } // namespace json

} // namespace detray