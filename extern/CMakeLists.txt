SET(ALGEBRA_PLUGIN_USE_VECMEM ON CACHE BOOL "Use Vecmem for array algebra")
add_subdirectory(algebra-plugins)

if (DETRAY_UNIT_TESTS)
    add_subdirectory(googletest)
endif()

if (DETRAY_IO_CSV)
    add_subdirectory(dfelibs)
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()
