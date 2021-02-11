/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "masks/mask_identifier.hpp"
#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "tools/planar_intersector.hpp"

#include <cmath>
#include <climits>

namespace detray
{
    /** This is a 2-dimensional mask for the annulus geometry that is
     *  e.g. used for the itk strip endcaps.
     * 
     * It is defined by the two radii _values[0] and  _values[1] in the polar
     * coordinate system of and endcap strip module, as well as the two phi 
     * boundaries, _values[2] and _values[3], that are only easily definable in 
     * the local strips system (relative to the average Phi along the 
     * strips r coordinate). Using a conversion between the two coordinate 
     * systems, these boundaries can be checked with a tolerance in r (t0 and 
     * t1), as well as phi (t3 and t3). 
     * Due to the local polar coordinate system of the strips needing a 
     * different origin from the discs polar one, three additional conversion
     * parameters are included (_values[4], values[5], _values[6]).
     * Here, the first two are the origin shift in xy, while _values[6] is the 
     * average Phi angle mentioned above.
     **/
    template <typename scalar_type,
              typename intersector_type = planar_intersector,
              typename links_type = bool,
              unsigned int kMaskIdentifier = e_annulus2>
    struct annulus2
    {
        using mask_values = darray<scalar_type, 7>;

        mask_values _values = {0., std::numeric_limits<scalar_type>::infinity(),
                               -std::numeric_limits<scalar_type>::infinity(), 
                               std::numeric_limits<scalar_type>::infinity(),
                               0., 0., 0.};

        links_type _links;

        static constexpr unsigned int mask_identifier = kMaskIdentifier;

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        annulus2<scalar_type, intersector_type, links_type, kMaskIdentifier> &
        operator=(const darray<scalar_type, 7> &rhs)
        {
            _values = rhs;
            return (*this);
        }

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked in local polar coord
         * @param t0 is the tolerance in minR
         * @param t1 is the tolerance in minPhi
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template<typename local_type>
        intersection_status is_inside(const typename local_type::point2 &p,
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
           // The two quantities to check: r^2 in module system, phi in strips system

           // In cartesian coordinates go to modules system by shifting origin
           if constexpr(std::is_same_v<local_type, __plugin::cartesian2>) {
              // Calculate radial coordinate in module system:
              scalar_type x_mod = p[0] - _values[4];
              scalar_type y_mod = p[1] - _values[5];
              scalar_type r_mod2 = x_mod * x_mod + y_mod * y_mod;

              // apply tolerances
              scalar_type minR_tol = _values[0] - t0;
              scalar_type maxR_tol = _values[1] + t0;

              if (r_mod2 < minR_tol*minR_tol or r_mod2 > maxR_tol*maxR_tol) return e_outside;

              scalar_type phi_strp = getter::phi(p) - _values[6];
              // Check phi boundaries, which are well def. in local frame
              return (phi_strp >= _values[2] - t1 and phi_strp <= _values[3] + t1) ? e_inside : e_outside;
            }
            // polar strip coordinates given
            else {
              // For a point p in local polar coordinates, rotate by avr phi
              scalar_type phi_strp = p[1] - _values[6];

              // Check phi boundaries, which are well def. in local frame
              if (phi_strp < _values[2] - t1 || phi_strp > _values[3] + t1) return e_outside;

              // Now go to module frame to check r boundaries. Use the origin shift
              // in polar coordinates for that
              typename local_type::point2 shift_xy = {-1*_values[4], -1*_values[5]};
              scalar_type shift_r   = getter::perp(shift_xy);
              scalar_type shift_phi = getter::phi(shift_xy);

              scalar_type r_mod2 = shift_r * shift_r + p[0] * p[0] + 
                                   2 * shift_r * p[0] * std::cos(phi_strp - shift_phi);

              // apply tolerances
              scalar_type minR_tol = _values[0] - t0;
              scalar_type maxR_tol = _values[1] + t0;

              return (r_mod2 >= minR_tol*minR_tol and r_mod2 <= maxR_tol*maxR_tol) ? e_inside : e_outside;
            }
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (_values == rhs);
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const annulus2<scalar_type> &rhs)
        {
            return operator==(rhs._values);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar_type &operator[](unsigned int value_index)
        {
            return _values[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar_type operator[](unsigned int value_index) const
        {
            return _values[value_index];
        }

        /** Return an associated intersector type */
        intersector_type intersector() const { return intersector_type{}; };

        /** Return the volume link - const reference */
        const links_type &links() const { return _links; }

        /** Return the volume link - non-const access */
        links_type &links() { return _links; }

    };

} // namespace detray
