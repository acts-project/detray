message(STATUS "Building 'detray::io' component")

if (DETRAY_IO_CSV)
    add_subdirectory(csv)
endif()
