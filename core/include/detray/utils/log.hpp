/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <iostream>

#ifdef DETRAY_ENABLE_LOGGING
#define DETRAY_DEBUG(x) std::cout << "DETRAY: " << x << std::endl;
#else
#define DETRAY_DEBUG(x)
#endif
