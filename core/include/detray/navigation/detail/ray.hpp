/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s).
#include <ostream>

namespace detray::detail {

/// @brief describes a straight-line trajectory
template <typename transform3_t>
class ray {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    using free_track_parameters_type = free_track_parameters<transform3_type>;
    using free_vector_type = typename free_track_parameters_type::vector_type;

    // Track helper
    using matrix_operator = typename transform3_t::matrix_actor;
    using track_helper = detail::track_helper<matrix_operator>;

    ray() = default;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    DETRAY_HOST_DEVICE ray(const point3 &pos, const scalar_type /*time*/,
                           const vector3 &dir, const scalar_type /*qop*/)
        : _pos{pos}, _dir{dir} {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    DETRAY_HOST_DEVICE ray(const free_vector_type &free_vec)
        : ray(track_helper().pos(free_vec), track_helper().time(free_vec),
              track_helper().dir(free_vec), track_helper().qop(free_vec)) {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    DETRAY_HOST_DEVICE ray(const free_track_parameters_type &track)
        : ray(track.vector()) {}

    /// @returns position on the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    /// @returns position on the ray paramterized by path length
    DETRAY_HOST_DEVICE point3 pos(const scalar_type s) const {
        // Direction is always normalized in the constructor
        return _pos + s * _dir;
    }

    /// @param position new position on the ray
    DETRAY_HOST_DEVICE void set_pos(point3 pos) { _pos = pos; }

    /// @returns direction of the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    /// @returns direction of the ray paramterized by path length
    DETRAY_HOST_DEVICE vector3 dir(const scalar_type /*s*/) const {
        return this->dir();
    }

    /// @param dir new direction of the ray
    DETRAY_HOST_DEVICE void set_dir(vector3 dir) { _dir = dir; }

    /// Print
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const ray &r) {
        os << "ray: ";
        os << "ori = [" << r._pos[0] << ", " << r._pos[1] << ", " << r._pos[2]
           << "], ";
        os << "dir = [" << r._dir[0] << ", " << r._dir[1] << ", " << r._dir[2]
           << "]" << std::endl;

        return os;
    }

    private:
    /// origin of ray
    point3 _pos{0.f, 0.f, 0.f};
    /// direction of ray
    vector3 _dir{0.f, 0.f, 1.f};
};

}  // namespace detray::detail
