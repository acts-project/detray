# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Helper function for setting up the detray libraries.
#
# Usage: detray_add_library( detray_core core "header1.hpp"... )
#
function( detray_add_library fullname basename )

   # Create the library.
   add_library( ${fullname} INTERFACE ${ARG_UNPARSED_ARGUMENTS} )

   # Set up how clients should find its headers.
   target_include_directories( ${fullname} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )

   # Make sure that the library is available as "detray::${basename}" in every
   # situation.
   set_target_properties( ${fullname} PROPERTIES EXPORT_NAME ${basename} )
   add_library( detray::${basename} ALIAS ${fullname} )

   # Set up the installation of the library and its headers.
   install( TARGETS ${fullname}
      EXPORT detray-exports
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" )
   install( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" )

endfunction( detray_add_library )

# Helper function for setting up the detray tests.
#
# Usage: detray_add_test( source1.cpp source2.cpp
#                         LINK_LIBRARIES detray::core )
#
function( detray_add_test name )

   # Parse the function's options.
   cmake_parse_arguments( ARG "" "" "LINK_LIBRARIES" ${ARGN} )

   # Create the test executable.
   set( test_exe_name "detray_test_${name}" )
   add_executable( ${test_exe_name} ${ARG_UNPARSED_ARGUMENTS} )
   if( ARG_LINK_LIBRARIES )
      target_link_libraries( ${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES} )
   endif()

   # Run the executable as the test.
   add_test( NAME ${test_exe_name}
      COMMAND ${test_exe_name} )

endfunction( detray_add_test )

# Helper function for adding individual flags to "flag variables".
#
# Usage: detray_add_flag( CMAKE_CXX_FLAGS "-Wall" )
#
function( detray_add_flag name value )

   # Escape special characters in the value:
   set( matchedValue "${value}" )
   foreach( c "*" "." "^" "$" "+" "?" )
      string( REPLACE "${c}" "\\${c}" matchedValue "${matchedValue}" )
   endforeach()

   # Check if the variable already has this value in it:
   if( "${${name}}" MATCHES "${matchedValue}" )
      return()
   endif()

   # If not, then let's add it now:
   set( ${name} "${${name}} ${value}" PARENT_SCOPE )

endfunction( detray_add_flag )
