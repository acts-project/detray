# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable( eigen_find_volume "eigen_find_volume.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::eigen )

detray_add_executable( eigen_intersect_all "eigen_intersect_all.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::eigen )

detray_add_executable( eigen_intersect_surfaces "eigen_intersect_surfaces.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::eigen )

detray_add_executable( eigen_masks "eigen_masks.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::eigen )
