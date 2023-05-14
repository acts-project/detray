# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the tutorial(s) that are algebra agnostic.
detray_add_tutorial( define_an_actor
   "define_an_actor.cpp"
   LINK_LIBRARIES detray::core )

# Set up the rest of the tutorials
set(source_files 
    "build_predefined_detectors.cpp"
    "detector_ray_scan.cpp"
    "do_it_yourself_detector.cpp"
    "navigation_inspection.cpp"
    "read_detector_from_file.cpp"
    "track_geometry_changes.cpp"
    "write_detector_to_file.cpp" )

foreach(file ${source_files})
   # Extract name of the tutorial
   string(REGEX REPLACE "\\.[^.]*$" "" file_stem ${file})

   detray_add_tutorial( ${file_stem}
      ${file}
      LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io 
                     detray::tutorial detray::test detray_tests_common detray::utils)
endforeach()
