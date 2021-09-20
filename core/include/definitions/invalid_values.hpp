/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

namespace detray {
    
    template< typename T >
    T invalid_value(){
	return T::invalid_value();
    }

    // specialization for int
    template<>
    int invalid_value(){
	return std::numeric_limits<int>::max();
    }

    // specialization for unsigned int    
    template<>
    unsigned int invalid_value(){
	return std::numeric_limits<unsigned int>::max();
    }

    // specialization for long unsigned int    
    template<>
    long unsigned int invalid_value(){
	return std::numeric_limits<long unsigned int>::max();
    }
    
} // namespace detray
