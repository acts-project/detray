/**
 * DETRAY library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Make sure the logging macros are available
#include "detray/utils/logging.hpp"

// Allow to temporarily disable macro based logging.
// Note: Has to be followed up by an include of "quiet_log_end"!
#if defined(__GNUC__) && !defined(__IN_QUIET_LOG_SECTION__)

#define __IN_QUIET_LOG_SECTION__

#pragma push_macro("DETRAY_VERBOSE_HOST")
#pragma push_macro("DETRAY_VERBOSE_DEVICE")
#pragma push_macro("DETRAY_VERBOSE_HOST_DEVICE")
#pragma push_macro("DETRAY_DEBUG_HOST")
#pragma push_macro("DETRAY_DEBUG_DEVICE")
#pragma push_macro("DETRAY_DEBUG_HOST_DEVICE")

#undef DETRAY_VERBOSE_HOST
#undef DETRAY_VERBOSE_DEVICE
#undef DETRAY_VERBOSE_HOST_DEVICE
#undef DETRAY_DEBUG_HOST
#undef DETRAY_DEBUG_DEVICE
#undef DETRAY_DEBUG_HOST_DEVICE

#define DETRAY_VERBOSE_HOST(x)
#define DETRAY_VERBOSE_DEVICE(x)
#define DETRAY_VERBOSE_HOST_DEVICE(x)
#define DETRAY_DEBUG_HOST(x)
#define DETRAY_DEBUG_DEVICE(x)
#define DETRAY_DEBUG_HOST_DEVICE(x)

#endif
