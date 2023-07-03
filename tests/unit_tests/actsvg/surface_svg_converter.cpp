// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_writer.hpp"


// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/masks/cylinder2D.hpp"
#include "detray/masks/masks.hpp"
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/tracks.hpp"
#include <type_traits>

using namespace actsvg;

template <typename algebra_t = __plugin::transform3<detray::scalar>>
struct to_svg {

    using transform3 = algebra_t;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    // Visitor to the detector mask store that is called on the mask
    // collection that contains the mask (shape) type of the surface
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline svg::object operator()(
        const mask_group_t& mask_group, 
        const index_t& index) const {
        
        // Get the concrete mask instance for this surface
        const auto& m = mask_group[index];

        int red = 200;
        int green = 150;
        int blue = 50;
    
        views::x_y x_y_view;

        style::fill fill_style{{{red, green, blue}}};
        style::stroke stroke_style{{{red, green, blue}}};

        fill_style._fc._highlight = {"mouseover", "mouseout"};
        fill_style._fc._hl_rgb = {0, 255, 0};
        fill_style._fc._opacity = 0.5;

        // svg object

        // ring2D (without hole)
        if constexpr (std::is_same_v<decltype(m), const detray::mask<detray::ring2D<>>&>){
            return draw::cir("t0", {0., 0.}, m[m.get_shape().e_outer_r], fill_style, stroke_style);
        }

        //display::surface

        return draw::text("t0", {0., 0.}, {"Error"});
        
    }
};


int main(int argc, char* argv[]) {
    // create mask
    //constexpr detray::scalar r{100.f * detray::unit<detray::scalar>::mm};
    //constexpr detray::scalar hz{4.f * detray::unit<detray::scalar>::mm};
    //detray::mask<detray::cylinder2D<>> c{0u, r, -hz, hz};

    //auto polygon_2d = x_y_view(polygon);

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    vecmem::host_memory_resource host_mr;
    const toy_detector_t det = detray::create_toy_geometry(host_mr, 4, 3);
    // The 13th surface in the detector should be a disc
    const auto disc_descr = det.surface_lookup()[13];
    const auto disc_surface = detray::surface{det, disc_descr};
    const auto svg = disc_surface.visit_mask<to_svg<>>();

    style::stroke stroke_black = style::stroke();
    auto x_y_a = draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    svg::file file;
    file.add_object(svg);
    file.add_object(x_y_a);

    std::ofstream stream;
    stream.open("detray_actsvg.svg");
    stream << file;
    stream.close();
}