// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"

// System include(s)
#include <concepts>
#include <string_view>

namespace detray::io::concepts {

/// Concept for detray io reader backends
template <typename D, typename R>
concept reader_backend =
    requires(const R rb, detector_builder<typename D::metadata> det_builder,
             typename D::name_map names) {

    typename R::payload_type;

    { R::tag }
    ->std::same_as<const std::string_view&>;

    {
        R::template from_payload<D>(det_builder, names,
                                    typename R::payload_type())
    }
    ->std::same_as<void>;
};

/// Concept for detray io writer backends
template <typename D, typename W>
concept writer_backend = requires(const W wb, D det,
                                  typename D::name_map names) {

    typename W::payload_type;

    { W::tag }
    ->std::same_as<const std::string_view&>;

    { W::to_payload(det, names) }
    ->std::same_as<typename W::payload_type>;
};

}  // namespace detray::io::concepts
