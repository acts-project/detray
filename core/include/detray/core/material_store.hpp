/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/** A compile-time material list that provides the correct material description
 * to surface */
template <template <typename...> class tuple_t = dtuple,
          typename ID = unsigned int, typename... materials_t>
struct material_store {

    public:
    template <typename... Args>
    using tuple_type = tuple_t<Args...>;

    using material_tuple = tuple_type<materials_t...>;

    DETRAY_HOST
};

}  // namespace detray