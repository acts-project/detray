/**
 * DETRAY library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/utils/logging.hpp"

// Allow to temporarily disable macro based logging.
// Note: Has to be preceded by an include of "quiet_log_start"!
#if defined(__GNUC__) && defined(__IN_QUIET_LOG_SECTION__)

#undef DETRAY_VERBOSE_HOST
#undef DETRAY_VERBOSE_DEVICE
#undef DETRAY_VERBOSE_HOST_DEVICE
#undef DETRAY_DEBUG_HOST
#undef DETRAY_DEBUG_DEVICE
#undef DETRAY_DEBUG_HOST_DEVICE

#pragma pop_macro("DETRAY_VERBOSE_HOST")
#pragma pop_macro("DETRAY_VERBOSE_DEVICE")
#pragma pop_macro("DETRAY_VERBOSE_HOST_DEVICE")
#pragma pop_macro("DETRAY_DEBUG_HOST")
#pragma pop_macro("DETRAY_DEBUG_DEVICE")
#pragma pop_macro("DETRAY_DEBUG_HOST_DEVICE")

#undef __IN_QUIET_LOG_SECTION__

#endif
