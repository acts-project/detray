/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/qualifiers.hpp"

// detray utils
#include "detray/definitions/indexing.hpp"

namespace detray {

template <typename scalar_t, typename context_t = dindex>
class constant_magnetic_field {

    public:
    using vector3 = __plugin::vector3<scalar_t>;
    using point3 = __plugin::point3<scalar_t>;
    using context_type = context_t;

    DETRAY_HOST_DEVICE
    constant_magnetic_field(vector3 field) : _field(field) {}

    DETRAY_HOST_DEVICE
    const vector3 &get_field(point3 /*pos*/, context_t /*ctx*/) const {
        return _field;
    }

    vector3 _field;
};

}  // namespace detray