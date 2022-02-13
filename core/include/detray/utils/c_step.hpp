
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

namespace detray {

/** Constrained step struct, it is used to regulate the stepping
 *
 */
template <size_t kDIM, template <typename, std::size_t> class array_t = darray>
struct c_step {

    /** Possible step constraints */
    enum type : int {
        e_accuracy = 0,
        e_navigation = 1,
        e_user = 2,
        e_limit = 3
    };

    array_t<scalar, kDIM> _values = {std::numeric_limits<scalar>::max(), 0.,
                                     std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::max()};

    const scalar operator()() const {
        return std::min_element(_values.begin(), _values.end());
    }
};

}  // namespace detray
