# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable( vc_array_find_volume "vc_array_find_volume.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::vc_array )

detray_add_executable( vc_array_intersect_all "vc_array_intersect_all.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::vc_array )

detray_add_executable( vc_array_intersect_surfaces
   "vc_array_intersect_surfaces.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::vc_array )

detray_add_executable( vc_array_masks "vc_array_masks.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::vc_array )
