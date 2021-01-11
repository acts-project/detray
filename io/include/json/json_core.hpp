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

        /** Adding a volume to the json output format
         * 
         * @tparam volume_type is the type of the volume object and @param v the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename volume_type>
        njson write_volume(const volume_type &v)
        {
            njson vj;
            vj["name"] = v.name;
            return vj;
        }

        /** Adding a transform to the json output format
         * 
         * @tparam transform_type is the type of the transform object and @param t the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename transform_type>
        njson write_transform(const transform_type &t)
        {
            njson tj;
            const auto trl = t.translation();
            tj["translation"] = {trl[0], trl[1], trl[2]};

            const auto rot = t.rotation();
#ifndef __plugin_without_matrix_element_accessor
            tj["rotation"] = {rot(0, 0), rot(0, 1), rot(0, 2), rot(1, 0), rot(1, 1), rot(1, 2), rot(2, 0), rot(2, 1), rot(2, 2)};
#else
            tj["rotation"] = {rot[0][0], rot[0][1], rot[0][2], rot[1][0], rot[1][1], rot[1][2], rot[2][0], rot[2][1], rot[2][2]};
#endif

            return tj;
        }

        /** Adding a surface to the json output format
         * 
         * @tparam surface_type is the type of the surface object and @param s the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename surface_type>
        njson write_surface(const surface_type &s)
        {
            njson sj;
            sj["transform_index"] = s.transform();
            sj["mask_index"] = s.mask();
            sj["volume_index"] = s.volume();
            return sj;
        }

        /** Adding a detector to the json output format
         * 
         * @tparam detector_type is the type of the detector object and @param d the object
         * 
         * @return a valid json object with substructure
         **/
        template <typename detector_type>
        njson write_detector(const detector_type &d)
        {

            // Detector json --------
            njson dj;
            dj["name"] = d.name();

            const auto &volumes = d.volumes();

            // Volumne section --------
            njson vj;
            std::for_each(volumes.begin(), volumes.end(), [&](const auto &v) { vj.push_back(write_volume(v)); });
            dj["volumes"] = vj;

            // Portal section --------
            const auto &portal_transforms = d.portal_transforms();
            njson ptj;
            std::for_each(portal_transforms.begin(), portal_transforms.end(), [&](const auto &pt) { ptj.push_back(write_transform(pt)); });
            dj["portal_transforms"] = ptj;

            const auto &portal_surfaces = d.portal_surfaces();
            njson pj;
            std::for_each(portal_surfaces.begin(), portal_surfaces.end(), [&](const auto &p) { pj.push_back(write_surface(p)); });
            dj["surfaces"] = pj;

            const auto &portal_masks = d.portal_masks();
            dj["portal_masks"] = write_all_masks(portal_masks);

            // Surface section --------
            const auto &surface_transforms = d.surface_transforms();
            njson stj;
            std::for_each(surface_transforms.begin(), surface_transforms.end(), [&](const auto &st) { stj.push_back(write_transform(st)); });
            dj["surface_transforms"] = stj;

            const auto &surfaces = d.surfaces();
            njson sj;
            std::for_each(surfaces.begin(), surfaces.end(), [&](const auto &s) { sj.push_back(write_surface(s)); });
            dj["surfaces"] = sj;

            const auto &surface_masks = d.surface_masks();
            dj["surface_masks"] = write_all_masks(surface_masks);

            return dj;
        }

    } // namespace json

} // namespace detray
