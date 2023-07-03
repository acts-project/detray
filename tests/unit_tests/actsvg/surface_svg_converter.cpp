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


using namespace actsvg;

int main(int argc, char* argv[]) {
    
    int red = 200;
    int green = 150;
    int blue = 50;
    
    views::x_y x_y_view;

    style::fill fill_style{{{red, green, blue}}};
    style::stroke stroke_style{{{red, green, blue}}};
    style::stroke stroke_black = style::stroke();

    fill_style._fc._highlight = {"mouseover", "mouseout"};
    fill_style._fc._hl_rgb = {0, 255, 0};
    fill_style._fc._opacity = 0.5;

    auto x_y_a =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    constexpr detray::scalar r{100.f * detray::unit<detray::scalar>::mm};
    constexpr detray::scalar hz{4.f * detray::unit<detray::scalar>::mm};

    // create mask
    detray::mask<detray::cylinder2D<>> c{0u, r, -hz, hz};
    std::cout << c[c.get_shape().e_r];

    //auto polygon_2d = x_y_view(polygon);

    // svg object
    auto svg =
        draw::circle("t0", {0., 0.}, c[c.get_shape().e_r], fill_style, stroke_style);

    svg::file file;
    file.add_object(svg);
    file.add_object(x_y_a);

    std::ofstream stream;
    stream.open("detray_actsvg.svg");
    stream << file;
    stream.close();
}
