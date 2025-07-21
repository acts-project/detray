/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <iostream>

#ifdef DETRAY_ENABLE_LOGGING
#define __FILENAME__ \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
// #define DETRAY_DEBUG(lvl, x) std::cout << lvl<< ": " << x << std::endl;
#define DETRAY_LOG(lvl, x)                                                  \
    std::cout << __FILENAME__ << ":" << __LINE__ << " " << lvl << ": " << x \
              << std::endl;
#else
#define DETRAY_DEBUG(x)
#endif

#define DETRAY_DEBUG(x) DETRAY_LOG("DEBUG", x)
#define DETRAY_ERROR(x) DETRAY_LOG("ERROR", x)
