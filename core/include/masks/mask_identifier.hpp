/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <limits>

namespace detray {

	/** Enumerate different mask types for convenience
	**/
	enum mask_context : unsigned int {
		e_known_types = 6,
		e_rectangle2 = 0,
		e_trapezoid2 = 1,
		e_ring2 = 2,
		e_cylinder3 = 3,
		e_single3 = 4,
		e_annulus2 = 5,
		e_unknown = std::numeric_limits<unsigned int>::max(),
	};

} // namespace detray
