# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::io' component")

# Create the "main I/O library".
detray_add_library( detray_io io )

# Include all active components.
if( DETRAY_IO_CSV )
   add_subdirectory( csv )
endif()

if ( DETRAY_IO_JSON )
   add_subdirectory( json )
endif()