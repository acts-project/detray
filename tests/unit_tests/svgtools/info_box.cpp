/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/masks/masks.hpp"

#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    std::string title = "info box";
    actsvg::style::fill title_blue;
    title_blue._fc._rgb = {0, 0, 150};
    title_blue._fc._opacity = 0.8;
    actsvg::style::font title_font;
    title_font._size = 24;
    title_font._fc = actsvg::style::color{{255, 255, 255}};

    std::vector<std::string> info = {"this is an info box", "a = 101",
                                     "this line should define the length"};

    actsvg::style::fill info_blue;
    info_blue._fc._rgb = {0, 0, 150};
    info_blue._fc._opacity = 0.4;
    actsvg::style::font info_font;
    info_font._size = 18;

    actsvg::style::stroke stroke;

    auto p = actsvg::draw::circle("test_circle", {-100, 100}, 15., title_blue, stroke);

    auto t =
        actsvg::draw::text("test_instruction", {-100, 120},
                   {"Move the mouse cursor over the the circle"}, info_font);

    auto info_box = actsvg::draw::connected_info_box("test_box_", {100, 100}, title,
                                             title_blue, title_font, info,
                                             info_blue, info_font, stroke, p);

    //auto info_box = actsvg::draw::connected_text("ok", {100,100}, info, info_font, actsvg::style::transform{}, p);

    actsvg::svg::file ifile;
    ifile.add_object(p);
    ifile.add_object(t);
    ifile.add_object(info_box);

    std::ofstream tstream;
    tstream.open("test_core_infobox.svg");
    tstream << ifile;
    tstream.close();
}