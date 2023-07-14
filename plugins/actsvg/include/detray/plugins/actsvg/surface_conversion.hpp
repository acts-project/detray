#include <math.h>

#include <type_traits>
#include <vector>

// Project include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/mask_conversion.hpp"
#include "detray/plugins/actsvg/transform_conversion.hpp"

namespace detray::actsvg_visualization {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;

/// @brief A functor to get the proto surface.
///
/// @note Does not account for the surface surface transform.
///
/// @returns An actsvg proto surface of the mask.
struct to_proto_surface {

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline actsvg::proto::surface<point3_container> operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return convert_mask(m);
    }
};

/// @brief Calculates the proto surface of a surface.
///
/// @note Accounts for the transform of the surface.
///
/// @param d_surface The detray surface.
/// @param context The context.
///
/// @returns An actsvg proto surface representing the surface.
template <typename detector_t>
proto_surface convert_surface(
    const detray::surface<detector_t>& d_surface,
    const typename detector_t::geometry_context& context) {
    auto p_surface = d_surface.template visit_mask<to_proto_surface>();
    p_surface._transform = convert_transform(d_surface.transform(context));

    return p_surface;
}
}  // namespace detray::actsvg_visualization