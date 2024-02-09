/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
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
inline auto pick_random(container_t container) {
    const auto idx = static_cast<std::size_t>(std::rand()) % container.size();
    return container[idx];
}

// Black
inline constexpr std::array black{0, 0, 0};
inline constexpr std::array dark_grey{171, 171, 171};
inline constexpr std::array mortar{89, 89, 89};
inline constexpr std::array suva_grey{137, 137, 137};
inline constexpr std::array very_light_grey{207, 207, 207};

// Red tones
inline constexpr std::array red{255, 0, 0};
inline constexpr std::array cardinal{167, 51, 63};
inline constexpr std::array madder{165, 28, 48};
inline constexpr std::array auburn{167, 51, 63};
inline constexpr std::array burgundy{116, 18, 29};
inline constexpr std::array chocolate_cosmos{88, 12, 31};
inline constexpr std::array macaroni_and_cheese{255, 188, 121};
inline constexpr std::array pumpkin{255, 128, 14};
inline constexpr std::array tenne{200, 82, 0};

// Blue tones
inline constexpr std::array blue{0, 0, 255};
inline constexpr std::array celestial_blue{62, 146, 204};
inline constexpr std::array cerulean{0, 107, 164};
inline constexpr std::array lapis_lazulli{42, 98, 143};
inline constexpr std::array picton_blue{95, 158, 209};
inline constexpr std::array prussian_blue1{19, 41, 61};
inline constexpr std::array prussian_blue2{22, 50, 79};
inline constexpr std::array indigo_dye{24, 67, 90};
inline constexpr std::array sail{162, 200, 236};

// Green tones
inline constexpr std::array green{0, 255, 0};
inline constexpr std::array celadon{190, 230, 206};
inline constexpr std::array aquamarine1{188, 255, 219};
inline constexpr std::array aquamarine2{141, 255, 205};
inline constexpr std::array emerald{104, 216, 155};
inline constexpr std::array shamrock_green{79, 157, 105};

inline std::vector<actsvg::style::color> black_theme(
    const actsvg::scalar opacity) {
    return {{black, opacity}};
}

inline std::vector<actsvg::style::color> red_theme(
    const actsvg::scalar opacity) {
    return {{cardinal, opacity},
            {madder, opacity},
            {auburn, opacity},
            {burgundy, opacity},
            {chocolate_cosmos, opacity}};
}

inline std::vector<actsvg::style::color> blue_theme(
    const actsvg::scalar opacity) {
    return {{celestial_blue, opacity},
            {lapis_lazulli, opacity},
            {prussian_blue1, opacity},
            {prussian_blue2, opacity},
            {indigo_dye, opacity}};
}

inline std::vector<actsvg::style::color> green_theme(
    const actsvg::scalar opacity) {
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
