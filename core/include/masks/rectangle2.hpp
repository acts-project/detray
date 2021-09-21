/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "masks/mask_identifier.hpp"
#include "tools/planar_intersector.hpp"

#include <climits>
#include <cmath>
#include <sstream>
#include <string>

namespace detray {

/** This is a simple 2-dimensional mask for a regular rectangle
 *
 * @tparam intersector_type is a struct used for intersecting this cylinder
 * @tparam local_type is the default local frame definition type
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by half length in local0 coordinates _values[0] and _values[1],
 * and can be checked with a tolerance in t[0] and t[1].
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <typename intersector_type = planar_intersector,
          typename local_type = __plugin::cartesian2,
          typename links_type = bool, unsigned int kMaskContext = e_rectangle2,
          template <typename, unsigned int> class array_type = darray>
struct rectangle2 {
  using mask_tolerance = array_type<scalar, 2>;

  using mask_values = array_type<scalar, 2>;

  using mask_links_type = links_type;

  mask_values _values = {std::numeric_limits<scalar>::infinity(),
                         std::numeric_limits<scalar>::infinity()};

  links_type _links;

  static constexpr unsigned int mask_context = kMaskContext;

  static constexpr unsigned int mask_indentifier = e_rectangle2;

  static constexpr mask_tolerance within_epsilon = {
      std::numeric_limits<scalar>::epsilon(),
      std::numeric_limits<scalar>::epsilon()};

  /** Assignment operator from an array, convenience function
   *
   * @param rhs is the right hand side object
   **/
  rectangle2<intersector_type, local_type, links_type, kMaskContext> &
  operator=(const array_type<scalar, 2> &rhs) {
    _values = rhs;
    return (*this);
  }

  /** Mask operation
   *
   * @tparam inside_local_type is the local type for inside checking
   *
   * @param p the point to be checked
   * @param t is the tolerance tuple in (l0, l1)
   *
   * @return an intersection status e_inside / e_outside
   **/
  template <typename inside_local_type>
  intersection_status
  is_inside(const point2 &p, const mask_tolerance &t = within_epsilon) const {
    return (std::abs(p[0]) <= _values[0] + t[0] and
            std::abs(p[1]) <= _values[1] + t[1])
               ? e_inside
               : e_outside;
  }

  /** Equality operator from an array, convenience function
   *
   * @param rhs is the rectangle to be compared with
   *
   * checks identity within epsilon and @return s a boolean*
   **/
  bool operator==(const array_type<scalar, 2> &rhs) { return (_values == rhs); }

  /** Equality operator
   *
   * @param rhs is the rectangle to be compared with
   *
   * checks identity within epsilon and @return s a boolean*
   **/
  bool operator==(const rectangle2<> &rhs) { return operator==(rhs._values); }

  /** Access operator - non-const
   * @return the reference to the member variable
   */
  scalar &operator[](unsigned int value_index) { return _values[value_index]; }

  /** Access operator - non-const
   * @return a copy of the member variable
   */
  scalar operator[](unsigned int value_index) const {
    return _values[value_index];
  }

  /** Return the values */
  const mask_values &values() const { return _values; }

  /** Return an associated intersector type */
  intersector_type intersector() const { return intersector_type{}; };

  /** Return the local frame type */
  local_type local() const { return local_type{}; }

  /** Return the volume link - const reference */
  const links_type &links() const { return _links; }

  /** Return the volume link - non-const access */
  links_type &links() { return _links; }

  /** Transform to a string for output debugging */
  std::string to_string() const {
    std::stringstream ss;
    ss << "rectangle2," << kMaskContext;
    for (const auto &v : _values) {
      ss << "," << v;
    }
    return ss.str();
  }
};

} // namespace detray
