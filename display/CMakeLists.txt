# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::display' component")

# Set up the library.
detray_add_library( detray_display display
   "include/detray/style/styles.hpp"
   "include/detray/view/draw.hpp"
   "include/detray/view/views.hpp" )
target_link_libraries( detray_display
   INTERFACE detray::core Matplot++::matplot )
