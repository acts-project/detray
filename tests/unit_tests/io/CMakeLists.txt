# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Helper function for running a test in a specific directory.
function( _run_test_in_dir test dir )
   # Make the directory.
   file( MAKE_DIRECTORY "${dir}" )
   # Run the test in that directory.
   set_tests_properties( "detray_test_${test}" PROPERTIES
      WORKING_DIRECTORY "${dir}" )
endfunction()

# Set up the I/O tests. Running the reader and writer tests in separate
# directories, so that they would not interfere with each other.
detray_add_test( io
   "io_json_payload.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::io_array )

detray_add_test( io_writer
   "io_json_detector_writer.cpp"
   LINK_LIBRARIES GTest::gtest_main vecmem::core detray::core_array detray::io_array detray::utils_array )
_run_test_in_dir( io_writer
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/io_writer_test_rundir" )

detray_add_test( io_roundtrip
    "io_json_detector_roundtrip.cpp"
    LINK_LIBRARIES GTest::gtest_main vecmem::core detray::core_array
    detray::io_array detray::test detray_tests_common detray::utils_array)
_run_test_in_dir( io_roundtrip
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/io_roundtrip_test_rundir" )
