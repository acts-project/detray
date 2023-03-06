# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable( smatrix_find_volume "smatrix_find_volume.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::smatrix )

detray_add_executable( smatrix_intersect_all "smatrix_intersect_all.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::smatrix )

detray_add_executable( smatrix_intersect_surfaces
   "smatrix_intersect_surfaces.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::smatrix )

detray_add_executable( smatrix_masks "smatrix_masks.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::smatrix )
