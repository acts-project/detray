/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// GTest include(s)
#include <gtest/gtest.h>

namespace detray::detail {

template <template <typename> class check_t, typename detector_t,
          typename config_t = typename check_t<detector_t>::config>
void register_checks(const detector_t &det,
                     const typename detector_t::name_map &vol_names,
                     const config_t &cfg = {}) {
    ::testing::RegisterTest(
        "detray_validation", cfg.name().c_str(), nullptr, cfg.name().c_str(),
        __FILE__, __LINE__,
        [=]() -> typename check_t<detector_t>::fixture_type * {
            return new check_t<detector_t>(det, vol_names, cfg);
        });
}

}  // namespace detray::detail
