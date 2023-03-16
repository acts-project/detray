# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable( array_find_volume "array_find_volume.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::array )

detray_add_executable( array_intersect_all "array_intersect_all.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::array )

detray_add_executable( array_intersect_surfaces "array_intersect_surfaces.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::array )

detray_add_executable( array_masks "array_masks.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::array )

detray_add_executable( array_propagation "array_propagation.cpp"
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::array)
if( OpenMP_CXX_FOUND )
  target_link_libraries( detray_array_propagation PRIVATE OpenMP::OpenMP_CXX )
endif()
