# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s) that are algebra agnostic.
#detray_add_example( cpu
#   "tools_actor_chain.cpp"
#   LINK_LIBRARIES detray::core )


# Build the example executable using the array plugin.
detray_add_example( track_geometry_changes
   "track_geometry_changes.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( build_predefined_detectors
   "build_predefined_detectors.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( write_detector_to_file
   "write_detector_to_file.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( read_detector_from_file
   "read_detector_from_file.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( detector_ray_scan
   "detector_ray_scan.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( define_an_actor
   "define_an_actor.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
