add_subdirectory(algebra-plugins)

if (DETRAY_UNIT_TESTS AND NOT gtest)
    add_subdirectory(googletest)
endif()

if (DETRAY_IO_CSV AND NOT dfelibs)
    add_subdirectory(dfelibs)
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()
