if (DETRAY_UNIT_TESTS)
    add_subdirectory(googletest)
endif()

if (DETRAY_IO_CSV)
    add_subdirectory(dfelibs)
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()