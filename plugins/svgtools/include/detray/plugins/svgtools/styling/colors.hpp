/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Actsvg include(s)
#include "actsvg/core/style.hpp"

// System include(s)
#include <array>
#include <vector>

namespace detray::svgtools::styling::colors {

/// @brief Picks a random element in the container.
template <typename container_t>
auto pick_random(container_t container) {
    const auto idx = static_cast<std::size_t>(std::rand()) % container.size();
    return container[idx];
}

// Black
constexpr std::array black{0, 0, 0};
constexpr std::array dark_grey{171, 171, 171};
constexpr std::array mortar{89, 89, 89};
constexpr std::array suva_grey{137, 137, 137};
constexpr std::array very_light_grey{207, 207, 207};

// Red tones
constexpr std::array red{255, 0, 0};
constexpr std::array cardinal{167, 51, 63};
constexpr std::array madder{165, 28, 48};
constexpr std::array auburn{167, 51, 63};
constexpr std::array burgundy{116, 18, 29};
constexpr std::array chocolate_cosmos{88, 12, 31};
constexpr std::array macaroni_and_cheese{255, 188, 121};
constexpr std::array pumpkin{255, 128, 14};
constexpr std::array tenne{200, 82, 0};

// Blue tones
constexpr std::array blue{0, 0, 255};
constexpr std::array celestial_blue{62, 146, 204};
constexpr std::array cerulean{0, 107, 164};
constexpr std::array lapis_lazulli{42, 98, 143};
constexpr std::array picton_blue{95, 158, 209};
constexpr std::array prussian_blue1{19, 41, 61};
constexpr std::array prussian_blue2{22, 50, 79};
constexpr std::array indigo_dye{24, 67, 90};
constexpr std::array sail{162, 200, 236};

// Green tones
constexpr std::array green{0, 255, 0};
constexpr std::array celadon{190, 230, 206};
constexpr std::array aquamarine1{188, 255, 219};
constexpr std::array aquamarine2{141, 255, 205};
constexpr std::array emerald{104, 216, 155};
constexpr std::array shamrock_green{79, 157, 105};

std::vector<actsvg::style::color> black_theme(const actsvg::scalar opacity) {
    return {{black, opacity}};
}

std::vector<actsvg::style::color> red_theme(const actsvg::scalar opacity) {
    return {{cardinal, opacity},
            {madder, opacity},
            {auburn, opacity},
            {burgundy, opacity},
            {chocolate_cosmos, opacity}};
}

std::vector<actsvg::style::color> blue_theme(const actsvg::scalar opacity) {
    return {{celestial_blue, opacity},
            {lapis_lazulli, opacity},
            {prussian_blue1, opacity},
            {prussian_blue2, opacity},
            {indigo_dye, opacity}};
}

std::vector<actsvg::style::color> green_theme(const actsvg::scalar opacity) {
    return {{emerald, opacity},
            {shamrock_green, opacity},
            {celadon, opacity},
            {aquamarine1, opacity},
            {aquamarine2, opacity}};
}

// Same color circle that is used in matplot plugin
struct tableau_colorblind10 {
    static std::vector<actsvg::style::color> grey_tones(
        const actsvg::scalar opacity) {
        return {{dark_grey, opacity},
                {mortar, opacity},
                {suva_grey, opacity},
                {very_light_grey, opacity}};
    }
    static std::vector<actsvg::style::color> blue_tones(
        const actsvg::scalar opacity) {
        return {{cerulean, opacity}, {picton_blue, opacity}, {sail, opacity}};
    }
    static std::vector<actsvg::style::color> red_tones(
        const actsvg::scalar opacity) {
        return {{tenne, opacity},
                {pumpkin, opacity},
                {macaroni_and_cheese, opacity}};
    }
};

}  // namespace detray::svgtools::styling::colors
