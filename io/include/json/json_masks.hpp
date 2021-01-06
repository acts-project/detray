/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <algorithm>
#include <tuple>

#include <nlohmann/json.hpp>

using njson = nlohmann::json;

namespace detray
{

    namespace json
    {

        /** Adding a mask to the json output format
         * 
         * @tparam mask_type is the type of the mask object and @param m the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename mask_type>
        njson write_mask(const mask_type &m)
        {
            njson mj;
            mj["values"] = m._values;
            mj["links"] = m._links;
            return mj;
        }

        /** Loop over the masks of a specific type and write them
         * 
         * @tparam mask_container_type The type of the mask container
         * @param masks the container object
         * 
         * @return a valid json object with substructure
         **/
        template <typename mask_container_type>
        njson write_specific_masks(const mask_container_type &masks)
        {
            njson mj;
            if (not masks.empty())
            {
                mj["identifier"] = masks.begin()->mask_identifier;
                njson mjs;
                std::for_each(masks.begin(), masks.end(), [&](const auto &m) { mjs.push_back(write_mask(m)); });
                mj["instances"] = mjs;
            }
            return mj;
        }

        /** Write the typed mask container to json
         * 
         * @tparam mask_container_type the type of the tuple entries
         * @param masks the typed masked tuple
         * 
         * @return a valid json object with substructure
         **/
        template <typename... mask_container_type>
        njson write_all_masks(std::tuple<mask_container_type...> const &masks)
        {
            njson mj;
            std::apply(
                [&](mask_container_type const &...mask_args) {
                    ((mj.push_back(write_specific_masks(mask_args)), ...));
                },
                masks);
            return mj;
        }

    } // namespace json

} // namespace detray