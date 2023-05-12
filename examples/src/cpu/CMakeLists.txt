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
detray_add_example( geometry_linking_cpu
   "check_geometry_linking.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( predefined_detectors_cpu
   "predefined_detectors.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( detector_writer_cpu
   "detector_writer.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
detray_add_example( detector_reader_cpu
   "detector_reader.cpp"
   LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io_array 
                  detray::example detray::test detray_tests_common detray::utils)
