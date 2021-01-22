/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

namespace detray {

	/** Enumerate different mask types for convenience
	**/
	enum mask_identifier : unsigned int {
		e_rectangle2 = 0,
		e_trapezoid2 = 1,
		e_ring2 = 2,
		e_cylinder3 = 3,
		e_single3 = 4,
		e_annulus2 = 5,
	};

} // namespace detray
