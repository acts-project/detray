#pragma once
/*
// Project include(s)
#include "detray/plugins/actsvg_visualization/surface.hpp"
#include "detray/plugins/actsvg_visualization/volume.hpp"
#include "detray/plugins/actsvg_visualization/detector.hpp"

// System include(s)
#include <assert.h>
#include <vector>

namespace detray::actsvg_visualization {

const detector::detector_options default_options{
    .v_options = volume::volume_options{
        .s_options = surface::surface_options{
            {{{0, 128, 0}, 0.9}, {{34, 139, 34}, 0.9}, {{60, 179, 113}, 0.9}, {{85, 107, 47}, 0.9}, {{124, 252, 0}, 0.9}, {{154, 205, 50}, 0.9}}
        },
        .p_options = surface::portal_options{
            .s_options = surface::surface_options{
                {{{128, 0, 128}, 0.9}, {{138, 43, 226}, 0.9}, {{148, 0, 211}, 0.9}, {{160, 32, 240}, 0.9}, {{186, 85, 211}, 0.9}, {{218, 112, 214}, 0.9}}
            },
            .l_options = link::link_options{
                3.
            }
        }
    }

    
};

    class options_builder {
        public:

        options_builder() = delete;

        options_builder(detector::detector_options d_options = default_options)
        : _options{d_options} {}

        auto& get_options(){
            return _options;
        }

        private:
        detector::detector_options _options;
    };

}*/