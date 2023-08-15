#pragma once

// Project include(s)
#include "detray/masks/annulus2D.hpp"
#include "detray/tools/generators.hpp"

// System include(s)
#include <iostream>
#include <cmath>

namespace detray::svg::utils{

template <typename transform_t, typename mask_t> 
inline auto global_corners(const transform_t& transform, const mask_t& mask)
{
    auto ret = local_corners(mask);
    for (size_t i = 0; i < ret.size(); i++){
        ret[i] = transform.point_to_global(ret[i]);
    }
    return ret;
}

template <typename mask_t> 
inline auto local_corners(const mask_t& mask)
{
    return vertices<mask_t::point2, mask_t::point3>(mask, 0);
}

template <> 
inline auto local_corners<mask<annulus2D<>>>(const mask<annulus2D<>>& m)
{
    using point3_t = mask<annulus2D<>>::point3_t;
    using scalar_t = mask<annulus2D<>>::scalar_type;
    const auto corner_arr = m.get_shape().corners(m.values());
    dvector<point3_t> vertices;
    for (size_t i = 0; i < corner_arr.size(); i++){
        const scalar_t z{0};
        const scalar_t x{corner_arr[i*2] * math_ns::cos(corner_arr[i*2+1])};
        const scalar_t y{corner_arr[i*2] * math_ns::sin(corner_arr[i*2+1])};
        vertices.push_back(point3_t{x, y, z});
    }
    return vertices;
}

}