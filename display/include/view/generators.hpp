/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/masks.hpp"

namespace detray
{

    using transform3 = __plugin::transform3;
    using point3 = transform3::point3;

    /** Generate phi values
     *
     * @param start_phi is the start for the arc generation
     * @param end_phi is the end of the arc generation
     * @param lseg is the number of segments used to gnerate the arc
     * 
     * @return a vector of phi values for the arc
     */
    static inline dvector<scalar>
    phi_values(scalar start_phi, scalar end_phi, unsigned int lseg)
    {
        dvector<scalar> values;
        values.reserve(lseg + 1);
        scalar step_phi = (end_phi - start_phi) / lseg;
        for (unsigned int istep = 0; istep <= lseg; ++istep)
        {
            values.push_back(start_phi + istep * step_phi);
        }
        return values;
    }

    /** Generate vertices, specialized for masks: annulus2
     *
     * @note template types are simply forwarded to mask
     * 
     * @param sf is the surface that generates vertices
     * @param ls is the number of line segments if
     *
     * @return a generated list of vertices
     */
    template <typename intersector_type,
              typename local_type,
              typename links_type,
              unsigned int kMaskContext>
    dvector<point3> vertices(const annulus2<intersector_type, local_type, links_type, kMaskContext> &annulus_mask, unsigned int lseg)
    {
        using point2 = __plugin::cartesian2::point2;

        const auto &m_values = annulus_mask.values();

        scalar min_r = m_values[0];
        scalar max_r = m_values[1];
        scalar min_phi_rel = m_values[2];
        scalar max_phi_rel = m_values[3];
        scalar avg_phi = m_values[4];
        scalar origin_x = m_values[5];
        scalar origin_y = m_values[6];

        point2 origin_m = {origin_x, origin_y};

        /// Helper method: find inner outer radius at edges in STRIP PC
        auto circIx = [](scalar O_x, scalar O_y, scalar r, scalar phi) -> point2 {
            //                      _____________________________________________
            //                     /      2  2                    2    2  2    2
            //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
            // x = --------------------------------------------------------------
            //                                  2
            //                                 m  + 1
            //
            // y = m*x
            //
            scalar m = std::tan(phi);
            point2 dir = {std::cos(phi), std::sin(phi)};
            scalar x1 = (O_x + O_y * m -
                         std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                                   2 * O_x * O_y * m - std::pow(O_y, 2) +
                                   std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                        (std::pow(m, 2) + 1);
            scalar x2 = (O_x + O_y * m +
                         std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                                   2 * O_x * O_y * m - std::pow(O_y, 2) +
                                   std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                        (std::pow(m, 2) + 1);

            point2 v1 = {x1, m * x1};
            if (vector::dot(v1, dir) > 0)
                return v1;
            return {x2, m * x2};
        };

        // calculate corners in STRIP XY, keep them we need them for minDistance()
        point2 ul_xy = circIx(origin_x, origin_y, max_r, max_phi_rel);
        point2 ll_xy = circIx(origin_x, origin_y, min_r, max_phi_rel);
        point2 ur_xy = circIx(origin_x, origin_y, max_r, min_phi_rel);
        point2 lr_xy = circIx(origin_x, origin_y, min_r, min_phi_rel);

        auto inner_phi = phi_values(getter::phi(ll_xy - origin_m), getter::phi(lr_xy - origin_m), lseg);
        auto outer_phi = phi_values(getter::phi(ur_xy - origin_m), getter::phi(ul_xy - origin_m), lseg);

        dvector<point3> annulus_vertices;
        annulus_vertices.reserve(inner_phi.size() + outer_phi.size());
        for (auto iphi : inner_phi)
        {
            annulus_vertices.push_back(point3{min_r * std::cos(iphi) + origin_x, min_r * std::sin(iphi) + origin_y, 0.});
        }

        for (auto ophi : outer_phi)
        {
            annulus_vertices.push_back(point3{max_r * std::cos(ophi) + origin_x, max_r * std::sin(ophi) + origin_y, 0.});
        }

        return annulus_vertices;
    }

    /** Generate vertices, specialized for masks: rectangle2
     *
     * @note template types are simply forwarded to mask
     *
     * @param sf is the surface that generates vertices
     * @param ls is the number of line segments (ignored for rectangles)
     *
     * @return a generated list of vertices
     */
    template <typename intersector_type,
              typename local_type,
              typename links_type,
              unsigned int kMaskContext>
    dvector<point3> vertices(const rectangle2<intersector_type, local_type, links_type, kMaskContext> &rectangle_mask, unsigned int /*ignored*/)
    {
        const auto &m_values = rectangle_mask.values();
        // left hand lower corner
        point3 lh_lc = {-m_values[0], -m_values[1], 0.};
        // right hand lower corner
        point3 rh_lc = {m_values[0], -m_values[1], 0.};
        // right hand upper corner
        point3 rh_uc = {m_values[0], m_values[1], 0.};
        // left hand upper corner
        point3 lh_uc = {-m_values[0], m_values[1], 0.};
        return {lh_lc, rh_lc, rh_uc, lh_uc};
        // Return the confining vertices
    }

    /** Generate vertices, specialized for masks: ring2
     *
     * @note template types are simply forwarded to mask
     *
     * @param sf is the surface that generates vertices
     * @param ls is the number of line segments if
     *
     * @return a generated list of vertices
     */
    template <typename intersector_type,
              typename local_type,
              typename links_type,
              unsigned int kMaskContext>
    dvector<point3> vertices(const ring2<intersector_type, local_type, links_type, kMaskContext> &ring_mask, unsigned int lseg)
    {
        return {};
    }

    /** Generate vertices, specialized for masks: trapezoid2
     *
     * @note template types are simply forwarded to mask
     *
     * @param sf is the surface that generates vertices
     * @param ls is the number of line segments if, ignored
     * 
     * @return a generated list of vertices
     */
    template <typename intersector_type,
              typename local_type,
              typename links_type,
              unsigned int kMaskContext>
    dvector<point3> vertices(const trapezoid2<intersector_type, local_type, links_type, kMaskContext> &trapezoid_mask, unsigned int /* ignored */)
    {

        const auto &m_values = trapezoid_mask.values();
        // left hand lower corner
        point3 lh_lc = {-m_values[0], -m_values[2], 0.};
        // right hand lower corner
        point3 rh_lc = {m_values[0], -m_values[2], 0.};
        // right hand upper corner
        point3 rh_uc = {m_values[1], m_values[2], 0.};
        // left hand upper corner
        point3 lh_uc = {-m_values[1], m_values[2], 0.};
        // Return the confining vertices
        return {lh_lc, rh_lc, rh_uc, lh_uc};
    }

} // namespace detray
