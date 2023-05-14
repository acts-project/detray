# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the tutorial(s) that are algebra agnostic.
detray_add_tutorial( define_an_actor
   "propagation/define_an_actor.cpp"
   LINK_LIBRARIES detray::core )

# Set up the rest of the tutorials
set(source_files 
    "detector/build_predefined_detectors.cpp"
    "detector/detector_ray_scan.cpp"
    "detector/do_it_yourself_detector.cpp"
    "detector/read_detector_from_file.cpp"
    "detector/track_geometry_changes.cpp"
    "detector/write_detector_to_file.cpp"
    "propagation/navigation_inspection.cpp" )

foreach(file ${source_files})
   # Extract name of the tutorial
   #string(REGEX REPLACE "\\.[^.]*$" "" file_stem ${file})
   get_filename_component( file_name ${file} NAME_WE )

   detray_add_tutorial( ${file_name}
      ${file}
      LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                     detray::tutorial detray::test detray_tests_common detray::utils)
endforeach()
