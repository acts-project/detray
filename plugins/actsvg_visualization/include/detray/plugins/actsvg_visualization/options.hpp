#pragma once

// Project include(s)
#include "detray/plugins/actsvg_visualization/surface.hpp"
#include "detray/plugins/actsvg_visualization/volume.hpp"
#include "detray/plugins/actsvg_visualization/detector.hpp"

// System include(s)
#include <assert.h>
#include <vector>

namespace detray::actsvg_visualization {

constexpr detector::detector_options default_options{};

    class options_builder {
        public:

        options_builder() = delete;

        options_builder(detector::detector_options d_options = default_options)
        : _options{d_options} {}

        auto& get_options(){
            return _options;
        }

        /*auto& portal_visible(bool value){
            return *this;
        }

        auto& surface_visible(bool value){
            return *this;
        }

        auto& link_visible(bool value){
            return *this;
        }*/

        private:
        detector::detector_options _options;
    };

}