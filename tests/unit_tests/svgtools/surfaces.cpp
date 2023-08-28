/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    // Creating the views.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the svg generator for the detector.
    const detray::svgtools::illustrator il{det, context};

    // Indexes of the surfaces in the detector to be visualized.
    std::array indices{200UL, 201UL, 202UL, 203UL, 204UL, 205UL, 206UL, 207UL, 208UL, 209UL,
                                          210UL, 211UL, 212UL, 213UL, 214UL, 215UL, 216UL, 217UL, 218UL, 219UL,
                                          220UL, 221UL, 222UL, 223UL, 224UL, 225UL, 226UL, 227UL, 228UL, 229UL,
                                          230UL, 231UL, 232UL, 233UL, 234UL, 235UL, 236UL, 237UL, 238UL, 239UL,
                                          240UL, 241UL, 242UL, 243UL, 244UL, 245UL, 246UL, 247UL, 248UL, 249UL,
                                          250UL, 251UL, 252UL, 253UL, 254UL, 255UL, 256UL, 257UL, 258UL, 259UL,
                                          260UL, 261UL, 262UL, 263UL, 264UL, 265UL, 266UL, 267UL, 268UL, 269UL,
                                          270UL, 271UL, 272UL, 273UL, 274UL, 275UL, 276UL, 277UL, 278UL, 279UL,
                                          280UL, 281UL, 282UL, 283UL, 284UL, 285UL, 286UL, 287UL, 288UL, 289UL,
                                          290UL, 291UL, 292UL, 293UL, 294UL, 295UL, 296UL, 297UL, 298UL, 299UL,
                                          300UL, 301UL, 302UL, 303UL, 304UL, 305UL, 306UL, 307UL, 308UL, 309UL,
                                          310UL, 311UL, 312UL, 313UL, 314UL, 315UL, 316UL, 317UL, 318UL, 319UL,
                                          320UL, 321UL, 322UL, 323UL, 324UL, 325UL, 326UL, 327UL, 328UL, 329UL,
                                          330UL, 331UL, 332UL, 333UL, 334UL, 335UL, 336UL, 337UL, 338UL, 339UL,
                                          340UL, 341UL, 342UL, 343UL, 344UL, 345UL, 346UL, 347UL, 348UL, 349UL,
                                          350UL, 351UL, 352UL, 353UL, 354UL, 355UL, 356UL, 357UL, 358UL, 359UL,
                                          360UL, 361UL, 362UL, 363UL, 364UL, 365UL, 366UL, 367UL, 368UL, 369UL,
                                          370UL, 371UL, 372UL, 373UL, 374UL, 375UL, 376UL, 377UL, 378UL, 379UL,
                                          380UL, 381UL, 382UL, 383UL, 384UL, 385UL, 386UL, 387UL, 388UL, 389UL,
                                          390UL, 391UL, 392UL, 393UL, 394UL, 395UL, 396UL, 397UL, 398UL, 399UL,
                                          400UL};

    for (std::size_t i : indices) {
        std::string name = "test_svgtools_surface" + std::to_string(i);
        // Visualization of portal i:
        const auto svg_xy = il.draw_surface(name, i, xy);
        detray::svgtools::write_svg(name + "_xy.svg", {axes, svg_xy});
        const auto svg_zr = il.draw_surface(name, i, zr);
        detray::svgtools::write_svg(name + "_zr.svg", {axes, svg_zr});
    }
}