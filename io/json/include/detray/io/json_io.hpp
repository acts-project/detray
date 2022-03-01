/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// This header is used to disable GCC errors that could be occuring in
// nlohmann/json.hpp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#include <nlohmann/json.hpp>
#pragma GCC diagnostic pop

namespace ALGEBRA_TRANSFORM3_NAMESPACE {

void to_json(nlohmann::json& j, const __plugin::transform3<detray::scalar>& t) {

    // Get translation and write it
    auto translation = t.translation();
    detray::darray<detray::scalar, 3> tdata = {translation[0], translation[1],
                                               translation[2]};
    j["translation"] = tdata;
    // Get rotation and write it
    auto rotation = t.rotation();
    detray::darray<detray::scalar, 9> rdata = {
        algebra::getter::element(rotation, 0, 0),
        algebra::getter::element(rotation, 0, 1),
        algebra::getter::element(rotation, 0, 2),
        algebra::getter::element(rotation, 1, 0),
        algebra::getter::element(rotation, 1, 1),
        algebra::getter::element(rotation, 1, 2),
        algebra::getter::element(rotation, 2, 0),
        algebra::getter::element(rotation, 2, 1),
        algebra::getter::element(rotation, 2, 2)};
    j["rotation"] = rdata;
}

void from_json(const nlohmann::json& j, __plugin::transform3<detray::scalar>& t) {

    detray::darray<detray::scalar, 9> rdata = {1., 0., 0., 0., 1.,
                                               0., 0., 0., 1.};
    if (not j["rotation"].empty()) {
        rdata = j["rotation"];
    }
    detray::darray<detray::scalar, 3> tdata = {0., 0., 0.};
    if (not j["translation"].empty()) {
        tdata = j["translation"];
    }

    detray::matrix<detray::scalar, 4, 4> ma;
    algebra::getter::element(ma, 0, 0) = rdata[0];
    algebra::getter::element(ma, 0, 1) = rdata[1];
    algebra::getter::element(ma, 0, 2) = rdata[2];
    algebra::getter::element(ma, 0, 3) = tdata[0];
    algebra::getter::element(ma, 1, 0) = rdata[3];
    algebra::getter::element(ma, 1, 1) = rdata[4];
    algebra::getter::element(ma, 1, 2) = rdata[5];
    algebra::getter::element(ma, 1, 3) = tdata[1];
    algebra::getter::element(ma, 2, 0) = rdata[6];
    algebra::getter::element(ma, 2, 1) = rdata[7];
    algebra::getter::element(ma, 2, 2) = rdata[8];
    algebra::getter::element(ma, 2, 3) = tdata[2];
    algebra::getter::element(ma, 3, 0) = 0.;
    algebra::getter::element(ma, 3, 1) = 0.;
    algebra::getter::element(ma, 3, 2) = 0.;
    algebra::getter::element(ma, 3, 3) = 1.;

    t = __plugin::transform3<detray::scalar>(ma);
}

}  // namespace __plugin::math
