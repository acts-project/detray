/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/surface.hpp"

// System include(s).
#include <vector>

namespace detray::detail {

/// @brief function that retrieves access to the portals of a cylindrical volume
///
/// The portals are returned as vectors in the order [inner, outer, lower,
/// upper]
template <typename detector_t>
auto get_cylinder_portals(const detector_volume<detector_t> &vol) {

    using scalar_t = typename detector_t::scalar_type;

    std::vector<const typename detector_t::volume_type &> inner_pt{},
        outer_pt{}, lower_pt{}, upper_pr{};

    std::map<const typename detector_t::surface_type &, scalar_t> radii{0.f};
    std::map<const typename detector_t::surface_type &, scalar_t> z_pos{0.f};

    // Loop over all portals
    for (const auto &pt_desc : vol.portals()) {
        auto pt = surface{vol.detector(), pt_desc};
        const std::string name = pt.shape_name();

        if (name == "cylinder2D" || name == "concentric_cylinder2D") {
            radii[pt_desc] = pt.boundary(cylinder2D::e_r);
            z_pos[pt_desc] = pt.boundary(cylinder2D::e_n_half_z);
            z_pos[pt_desc] = pt.boundary(cylinder2D::e_p_half_z);
        } else {
            radii[pt_desc] = pt.boundary(ring2D::e_inner_r);
            radii[pt_desc] = pt.boundary(ring2D::e_outer_r);
        }
    }
}

}  // namespace detray::detail
