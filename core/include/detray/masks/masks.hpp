/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/masks/annulus2D.hpp"
#include "detray/masks/cuboid3D.hpp"
#include "detray/masks/cylinder2D.hpp"
#include "detray/masks/cylinder3D.hpp"
#include "detray/masks/line.hpp"
#include "detray/masks/rectangle2D.hpp"
#include "detray/masks/ring2D.hpp"
#include "detray/masks/single3D.hpp"
#include "detray/masks/trapezoid2D.hpp"
#include "detray/masks/unmasked.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <sstream>
#include <vector>

namespace detray {

/// @brief Mask a region on a surface and link it to a volume.
///
/// The class uses a lightweight 'shape' that defines the local geometry of a
/// surface (local coordinates/axes, the intersection algorithm, as well as how
/// to check boundaries on that surface). A simple example for such a surface
/// shape implementation is @c rectangle2D .
/// This class can then be instantiated as a concrete mask that holds the
/// boundary values, as well as the link to a particular volume, which is
/// needed during geometry navigation.
///
/// @tparam shape_t underlying geometrical shape of the mask, rectangle etc
/// @tparam links_t the type of link into the volume container
///                 (e.g. single index vs range)
template <typename shape_t, typename links_t = dindex,
          typename algebra_t = __plugin::transform3<detray::scalar>,
          template <typename, std::size_t> class array_t = darray>
class mask {
    public:
    using links_type = links_t;
    using scalar_type = typename algebra_t::scalar_type;
    using shape = shape_t;
    using boundaries = typename shape::boundaries;
    using mask_values = array_t<scalar_type, boundaries::e_size>;
    using local_frame_type =
        typename shape::template local_frame_type<algebra_t>;
    using measurement_frame_type =
        typename shape::template measurement_frame_type<algebra_t>;
    // Linear algebra types
    using loc_point_t = typename shape::template loc_point_type<algebra_t>;
    using point3_t = typename algebra_t::point3;
    using matrix_operator = typename algebra_t::matrix_actor;
    using size_type = typename algebra_t::size_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename algebra_t::template matrix_type<ROWS, COLS>;

    /// Default constructor
    constexpr mask() = default;

    /// Constructor from single mask boundary values
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr mask(const links_type& link,
                                               Args&&... args)
        : _values({{std::forward<Args>(args)...}}), _volume_link(link) {}

    /// Constructor from mask boundary vector
    DETRAY_HOST mask(const std::vector<scalar_type>& values,
                     const links_type& link)
        : _volume_link(link) {
        assert(values.size() == boundaries::e_size &&
               " Given number of boundaries does not match mask shape.");
        std::copy(std::cbegin(values), std::cend(values), std::begin(_values));
    }

    /// Assignment operator from an array, convenience function
    ///
    /// @param rhs is the right hand side object
    DETRAY_HOST
    auto operator=(const mask_values& rhs)
        -> mask<shape_t, links_t, algebra_t, array_t>& {
        _values = rhs;
        return (*this);
    }

    /// Equality operator
    ///
    /// @param rhs is the mask to be compared
    ///
    /// checks identity within epsilon and @returns a boolean if the values and
    /// links are equal.
    DETRAY_HOST_DEVICE
    bool operator==(
        const mask<shape_t, links_t, algebra_t, array_t>& rhs) const {
        return (_values == rhs._values && _volume_link == rhs._volume_link);
    }

    /// Access operator - non-const
    ///
    /// @returns the reference to the member variable
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const std::size_t value_index) -> scalar_type& {
        return _values[value_index];
    }

    /// Access operator - const
    ///
    /// @returns a copy of the member variable
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const std::size_t value_index) const
        -> scalar_type {
        return _values[value_index];
    }

    /// @returns the boundary values
    DETRAY_HOST_DEVICE
    inline constexpr auto get_shape() const -> const shape& { return _shape; }

    /// @returns the functor that projects a global cartesian point onto
    /// the local geometric coordinate system.
    template <typename transform3_t>
    DETRAY_HOST_DEVICE inline auto to_local_frame(
        const transform3_t& trf, const point3_t& glob_p,
        const point3_t& glob_dir = {}) const -> loc_point_t {
        return local_frame_type{}.global_to_local(trf, glob_p, glob_dir);
    }

    /// @returns the intersection functor for the underlying surface geometry.
    DETRAY_HOST_DEVICE
    inline constexpr auto intersector() const ->
        typename shape::template intersector_type<algebra_t> {
        return {};
    }

    /// @brief Mask this shape onto a surface.
    ///
    /// @note the point is expected to be given in global cartesian coordinates
    /// by the caller. For the projection from global to local coordinates, the
    /// shape type is used.
    ///
    /// @param loc_p the point to be checked in the local polar focal system
    /// @param tol dynamic tolerance determined by caller
    ///
    /// @return an intersection status e_inside / e_outside
    DETRAY_HOST_DEVICE
    inline auto is_inside(
        const loc_point_t& loc_p,
        const scalar_type t = std::numeric_limits<scalar_type>::epsilon()) const
        -> intersection::status {

        return _shape.check_boundaries(_values, loc_p, t)
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /// @returns return local frame object
    DETRAY_HOST_DEVICE inline constexpr local_frame_type local_frame() const {
        return local_frame_type{};
    }

    /// @returns the boundary values
    DETRAY_HOST_DEVICE
    auto values() const -> const mask_values& { return _values; }

    /// @returns the volume link - const reference
    DETRAY_HOST_DEVICE
    auto volume_link() const -> const links_type& { return _volume_link; }

    /// @returns the volume link - non-const access
    DETRAY_HOST_DEVICE
    auto volume_link() -> links_type& { return _volume_link; }

    /// @returns a string representation of the mask
    DETRAY_HOST
    auto to_string() const -> std::string {
        std::stringstream ss;
        ss << shape::name;
        for (const auto& v : _values) {
            ss << ", " << v;
        }
        return ss.str();
    }

    template <size_type parameter_dim>
    DETRAY_HOST_DEVICE matrix_type<2, parameter_dim> projection_matrix() const {

        auto ret = matrix_operator().template zero<2, parameter_dim>();
        for (unsigned int i = 0u; i < shape::meas_dim; i++) {
            matrix_operator().element(ret, i, i) = 1.f;
        }

        return ret;
    }

    private:
    shape _shape;
    mask_values _values;
    links_type _volume_link;
};

}  // namespace detray