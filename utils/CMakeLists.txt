# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::utils' component")

# Set up the library.
detray_add_library( detray_utils utils 
    # detectors include(s)
    "include/detray/detectors/create_telescope_detector.hpp"
    "include/detray/detectors/create_toy_geometry.hpp"
    "include/detray/detectors/detector_metadata.hpp"
    # simulation include(s)
    "include/detray/simulation/event_generator/random_track_generator.hpp"
    "include/detray/simulation/event_generator/track_generators.hpp"
    "include/detray/simulation/event_generator/uniform_track_generator.hpp"
    "include/detray/simulation/event_writer.hpp"
    "include/detray/simulation/measurement_smearer.hpp"
    "include/detray/simulation/random_scatterer.hpp"
    "include/detray/simulation/simulator.hpp"
)

target_link_libraries( detray_utils INTERFACE detray_core)
