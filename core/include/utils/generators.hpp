/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "masks/masks.hpp"

namespace detray {

using transform3 = __plugin::transform3;
using point2 = __plugin::point2;

/** Generate phi values
 *
 * @param start_phi is the start for the arc generation
 * @param end_phi is the end of the arc generation
 * @param lseg is the number of segments used to gnerate the arc
 *
 * @return a vector of phi values for the arc
 */
static inline dvector<scalar> phi_values(scalar start_phi, scalar end_phi,
                                         unsigned int lseg) {
    dvector<scalar> values;
    values.reserve(lseg + 1);
    scalar step_phi = (end_phi - start_phi) / lseg;
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
template <typename intersector_type, typename local_type, typename links_type,
          unsigned int kMaskContext>
dvector<point3> vertices(const annulus2<intersector_type, local_type,
                                        links_type, kMaskContext> &annulus_mask,
                         unsigned int lseg) {

    const auto &m_values = annulus_mask.values();

    scalar min_r = m_values[0];
    scalar max_r = m_values[1];
    scalar min_phi_rel = m_values[2];
    scalar max_phi_rel = m_values[3];
    // scalar avg_phi = m_values[4];
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
        scalar x1 =
            (O_x + O_y * m -
             std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                       std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                       std::pow(r, 2))) /
            (std::pow(m, 2) + 1);
        scalar x2 =
            (O_x + O_y * m +
             std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                       std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                       std::pow(r, 2))) /
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

    auto inner_phi = phi_values(getter::phi(ll_xy - origin_m),
                                getter::phi(lr_xy - origin_m), lseg);
    auto outer_phi = phi_values(getter::phi(ur_xy - origin_m),
                                getter::phi(ul_xy - origin_m), lseg);

    dvector<point3> annulus_vertices;
    annulus_vertices.reserve(inner_phi.size() + outer_phi.size());
    for (auto iphi : inner_phi) {
        annulus_vertices.push_back(point3{min_r * std::cos(iphi) + origin_x,
                                          min_r * std::sin(iphi) + origin_y,
                                          0.});
    }

    for (auto ophi : outer_phi) {
        annulus_vertices.push_back(point3{max_r * std::cos(ophi) + origin_x,
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
template <bool kRadialCheck, typename intersector_type, typename local_type,
          typename links_type, unsigned int kMaskContext>
dvector<point3> vertices(
    const cylinder3<kRadialCheck, intersector_type, local_type, links_type,
                    kMaskContext> &annulus_mask,
    unsigned int lseg) {

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
template <typename intersector_type, typename local_type, typename links_type,
          unsigned int kMaskContext>
dvector<point3> vertices(
    const rectangle2<intersector_type, local_type, links_type, kMaskContext>
        &rectangle_mask,
    unsigned int /*ignored*/) {
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
template <typename intersector_type, typename local_type, typename links_type,
          unsigned int kMaskContext>
dvector<point3> vertices(const ring2<intersector_type, local_type, links_type,
                                     kMaskContext> &ring_mask,
                         unsigned int lseg) {
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
template <typename intersector_type, typename local_type, typename links_type,
          unsigned int kMaskContext>
dvector<point3> vertices(
    const trapezoid2<intersector_type, local_type, links_type, kMaskContext>
        &trapezoid_mask,
    unsigned int /* ignored */) {

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

/** Specialized method to generate vertices per maks group
 *
 * @tparam mask_group is the type of the split out mask group from all masks
 * @tparam mask_range is the type of the according mask range object
 *
 * @param masks is the associated (and split out) mask group
 * @param range is the range list of masks to be processed
 *
 * @return a jagged vector of points of the mask vertices (one per maks)
 **/
template <typename mask_group, typename mask_range>
auto vertices_for_mask_group(const mask_group &masks, const mask_range &range,
                             dindex n_segments = 1) {
    dvector<dvector<point3>> mask_vertices = {};
    for (auto i : sequence(range)) {
        const auto &mask = masks[i];
        mask_vertices.push_back(vertices(mask, n_segments));
    }
    return mask_vertices;
}

/** Variadic unrolled intersection - last entry
 *
 * @tparam mask_container is the type of the type of the mask container
 * @tparam mask_range is the mask range type
 * @tparam last_mask_context is the last mask group context
 *
 * @param masks the masks container
 * @param range the range within the mask group to be checked
 * @param mask_context the last mask group context
 *
 * @return a jagged vector of points of the mask vertices (one per maks)
 **/
template <typename mask_container, typename mask_range,
          dindex last_mask_context>
auto vertices_for_last_mask_group(const mask_container &masks,
                                  const mask_range &range,
                                  dindex mask_context) {
    dvector<dvector<point3>> mask_vertices = {};
    if (mask_context == last_mask_context) {
        mask_vertices = vertices_for_mask_group(
            masks.template group<last_mask_context>(), range);
    }
    return mask_vertices;
}

/** Variadic unrolled intersection - any integer sequence
 *
 * @tparam mask_container is the type of the type of the mask container
 * @tparam mask_range is the mask range type
 * @tparam first_mask_context is the first mask group context
 *
 * @param masks the masks container
 * @param range the range within the mask group to be checked
 * @param mask_context the last mask group context
 * @param available_contices the mask contices to be checked
 *
 * @return a jagged vector of points of the mask vertices (one per maks)
 **/
template <typename mask_container, typename mask_range,
          dindex first_mask_context, dindex... remaining_mask_context>
auto unroll_masks_for_vertices(
    const mask_container &masks, const mask_range &range, dindex mask_context,
    std::integer_sequence<dindex, first_mask_context, remaining_mask_context...>
        available_contices) {
    // Pick the first one for interseciton
    if (mask_context == first_mask_context) {
        return vertices_for_mask_group(
            masks.template group<first_mask_context>(), range);
    }
    // The reduced integer sequence
    std::integer_sequence<dindex, remaining_mask_context...> remaining;
    // Unroll as long as you have at least 2 entries
    if constexpr (remaining.size() > 1) {
        auto mask_vertices =
            unroll_masks_for_vertices(masks, range, mask_context, remaining);
        if (not mask_vertices.empty()) {
            return mask_vertices;
        }
    }
    // Last chance - intersect the last index if possible
    return vertices_for_last_mask_group<
        mask_container, mask_range,
        std::tuple_size_v<typename mask_container::mask_tuple> - 1>(
        masks, range, mask_context);
}

/** Create a r-phi polygon from principle parameters
 *
 * @param rmin minum r parameter
 * @param rmax maximum r parameter
 * @param phimin minimum phi parameter
 * @param phimax maximum phi parameters
 *
 * @return a polygon representation of the bin
 **/
std::vector<point2> r_phi_polygon(scalar rmin, scalar rmax, scalar phimin,
                                  scalar phimax, unsigned int nsegments = 1) {

    std::vector<point2> r_phi_poly;
    r_phi_poly.reserve(2 * nsegments + 2);

    scalar cos_min_phi = std::cos(phimin);
    scalar sin_min_phi = std::sin(phimin);
    scalar cos_max_phi = std::cos(phimax);
    scalar sin_max_phi = std::sin(phimax);

    // @TODO add phi generators
    r_phi_poly.push_back({rmin * cos_min_phi, rmin * sin_min_phi});
    r_phi_poly.push_back({rmin * cos_max_phi, rmin * sin_max_phi});
    r_phi_poly.push_back({rmax * cos_max_phi, rmax * sin_max_phi});
    r_phi_poly.push_back({rmax * cos_min_phi, rmax * sin_min_phi});

    return r_phi_poly;
};

}  // namespace detray
