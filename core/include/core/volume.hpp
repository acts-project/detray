/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

namespace detray
{
    /** Templated volume class
     * 
     * @tparam surface_container_type type of the internal surface container
     * @tparam source_type the type of the volume representation 
     */
    template <typename surface_container_type,
              typename portal_container_type = int, 
              typename source_type = int>
    class volume {

        volume(surface_container_type&& surfaces) 
         : _surfaces(std::move(surfaces));

        template <typename candidate_container_type>
        candidate_container_type surface_candidates() const {
            return _surfaces.get_candidates();
        }
    };
}