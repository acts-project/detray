/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

namespace detray {

    template < typename T >
    T invalid_value(){
	return T::invalid_value();
    }

    template<>
    unsigned int invalid_value(){
	return std::numeric_limits<unsigned int>::max();
    }

    template<>
    int invalid_value(){
	return std::numeric_limits<int>::max();
    }
    
    template<>
    unsigned long invalid_value(){
	return std::numeric_limits<unsigned long>::max();
    }
    
} // namespace detray
