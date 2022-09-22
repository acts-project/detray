/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/masks/masks.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

/** Generate phi values
 *
 * @param start_phi is the start for the arc generation
 * @param end_phi is the end of the arc generation
 * @param lseg is the number of segments used to gnerate the arc
 *
 * @return a vector of phi values for the arc
 */
template <typename scalar_t>
static inline dvector<scalar_t> phi_values(scalar_t start_phi, scalar_t end_phi,
                                           unsigned int lseg) {
    dvector<scalar_t> values;
    values.reserve(lseg + 1);
    scalar_t step_phi = (end_phi - start_phi) / lseg;
    for (unsigned int istep = 0; istep <= lseg; ++istep) {
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
template <typename point2_t, typename point3_t, typename links_t,
          typename transform3_t>
dvector<point3_t> vertices(
    const mask<annulus2D<>, links_t, transform3_t> &annulus_mask,
    unsigned int lseg) {
    using scalar_t = typename transform3_t::scalar_type;

    const auto &m_values = annulus_mask.values();

    scalar_t min_r = m_values[0];
    scalar_t max_r = m_values[1];
    scalar_t min_phi_rel = m_values[2];
    scalar_t max_phi_rel = m_values[3];
    // scalar_t avg_phi = m_values[4];
    scalar_t origin_x = m_values[5];
    scalar_t origin_y = m_values[6];

    point2_t origin_m = {origin_x, origin_y};

    /// Helper method: find inner outer radius at edges in STRIP PC
    auto circIx = [](scalar_t O_x, scalar_t O_y, scalar_t r,
                     scalar_t phi) -> point2_t {
        //                      _____________________________________________
        //                     /      2  2                    2    2  2    2
        //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
        // x = --------------------------------------------------------------
        //                                  2
        //                                 m  + 1
        //
        // y = m*x
        //
        scalar_t m = std::tan(phi);
        point2_t dir = {std::cos(phi), std::sin(phi)};
        scalar_t x1 =
            (O_x + O_y * m -
             std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                       std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                       std::pow(r, 2))) /
            (std::pow(m, 2) + 1);
        scalar_t x2 =
            (O_x + O_y * m +
             std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                       std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                       std::pow(r, 2))) /
            (std::pow(m, 2) + 1);

        point2_t v1 = {x1, m * x1};
        if (vector::dot(v1, dir) > 0)
            return v1;
        return {x2, m * x2};
    };

    // calculate corners in STRIP XY, keep them we need them for minDistance()
    point2_t ul_xy = circIx(origin_x, origin_y, max_r, max_phi_rel);
    point2_t ll_xy = circIx(origin_x, origin_y, min_r, max_phi_rel);
    point2_t ur_xy = circIx(origin_x, origin_y, max_r, min_phi_rel);
    point2_t lr_xy = circIx(origin_x, origin_y, min_r, min_phi_rel);

    auto inner_phi = phi_values(getter::phi(ll_xy - origin_m),
                                getter::phi(lr_xy - origin_m), lseg);
    auto outer_phi = phi_values(getter::phi(ur_xy - origin_m),
                                getter::phi(ul_xy - origin_m), lseg);

    dvector<point3_t> annulus_vertices;
    annulus_vertices.reserve(inner_phi.size() + outer_phi.size());
    for (auto iphi : inner_phi) {
        annulus_vertices.push_back(point3_t{min_r * std::cos(iphi) + origin_x,
                                            min_r * std::sin(iphi) + origin_y,
                                            0.});
    }

    for (auto ophi : outer_phi) {
        annulus_vertices.push_back(point3_t{max_r * std::cos(ophi) + origin_x,
                                            max_r * std::sin(ophi) + origin_y,
                                            0.});
    }

    return annulus_vertices;
}

/** Generate vertices, spacialized for masks: cylinder3
 *
 * @note template types are simply forwarded to mask
 *
 * @param sf is the surface that generates vertices
 * @param ls is the number of line segments if
 *
 * @return a generated list of vertices
 */
template <typename point2_t, typename point3_t, bool kRadialCheck,
          template <class> typename intersector_t, typename links_t,
          typename transform3_t>
dvector<point3_t> vertices(
    const mask<cylinder2D<kRadialCheck, intersector_t>, links_t, transform3_t>
        & /*cylinder_mask*/,
    unsigned int /*lseg*/) {

    return {};
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
template <typename point2_t, typename point3_t, typename links_t,
          typename transform3_t>
dvector<point3_t> vertices(
    const mask<rectangle2D<>, links_t, transform3_t> &rectangle_mask,
    unsigned int /*ignored*/) {
    const auto &m_values = rectangle_mask.values();
    // left hand lower corner
    point3_t lh_lc = {-m_values[0], -m_values[1], 0.};
    // right hand lower corner
    point3_t rh_lc = {m_values[0], -m_values[1], 0.};
    // right hand upper corner
    point3_t rh_uc = {m_values[0], m_values[1], 0.};
    // left hand upper corner
    point3_t lh_uc = {-m_values[0], m_values[1], 0.};
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
template <typename point2_t, typename point3_t, typename links_t,
          typename transform3_t>
dvector<point3_t> vertices(
    const mask<ring2D<>, links_t, transform3_t> & /*ring_mask*/,
    unsigned int /*lseg*/) {
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
template <typename point2_t, typename point3_t, typename links_t,
          typename transform3_t>
dvector<point3_t> vertices(
    const mask<trapezoid2D<>, links_t, transform3_t> &trapezoid_mask,
    unsigned int /* ignored */) {

    const auto &m_values = trapezoid_mask.values();
    // left hand lower corner
    point3_t lh_lc = {-m_values[0], -m_values[2], 0.};
    // right hand lower corner
    point3_t rh_lc = {m_values[0], -m_values[2], 0.};
    // right hand upper corner
    point3_t rh_uc = {m_values[1], m_values[2], 0.};
    // left hand upper corner
    point3_t lh_uc = {-m_values[1], m_values[2], 0.};
    // Return the confining vertices
    return {lh_lc, rh_lc, rh_uc, lh_uc};
}

/// Functor to produce vertices on a mask collection in a mask tuple container.
template <typename point2_t, typename point3_t>
struct vertexer {

    using output_type = dvector<dvector<point3_t>>;

    /// Specialized method to generate vertices per maks group
    ///
    /// @tparam mask_group_t is the type of the mask collection in a mask cont.
    /// @tparam mask_range_t is the type of the according mask range object
    ///
    /// @param masks is the associated (and split out) mask group
    /// @param range is the range list of masks to be processed
    ///
    /// @return a jagged vector of points of the mask vertices (one per maks)
    template <typename mask_group_t, typename mask_range_t>
    output_type operator()(const mask_group_t &masks, const mask_range_t &range,
                           dindex n_segments = 1) {
        output_type mask_vertices = {};
        for (auto i : sequence(range)) {
            const auto &mask = masks[i];
            mask_vertices.push_back(
                vertices<point2_t, point3_t>(mask, n_segments));
        }
        return mask_vertices;
    }
};

/** Create a r-phi polygon from principle parameters
 *
 * @param rmin minum r parameter
 * @param rmax maximum r parameter
 * @param phimin minimum phi parameter
 * @param phimax maximum phi parameters
 *
 * @return a polygon representation of the bin
 **/
template <typename scalar_t, typename point2_t>
std::vector<point2_t> r_phi_polygon(scalar_t rmin, scalar_t rmax,
                                    scalar_t phimin, scalar_t phimax,
                                    unsigned int nsegments = 1) {

    std::vector<point2_t> r_phi_poly;
    r_phi_poly.reserve(2 * nsegments + 2);

    scalar_t cos_min_phi = std::cos(phimin);
    scalar_t sin_min_phi = std::sin(phimin);
    scalar_t cos_max_phi = std::cos(phimax);
    scalar_t sin_max_phi = std::sin(phimax);

    // @TODO add phi generators
    r_phi_poly.push_back({rmin * cos_min_phi, rmin * sin_min_phi});
    r_phi_poly.push_back({rmin * cos_max_phi, rmin * sin_max_phi});
    r_phi_poly.push_back({rmax * cos_max_phi, rmax * sin_max_phi});
    r_phi_poly.push_back({rmax * cos_min_phi, rmax * sin_min_phi});

    return r_phi_poly;
}

}  // namespace detray