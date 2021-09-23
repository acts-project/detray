SET(ALGEBRA_PLUGIN_USE_VECMEM ON CACHE BOOL "Use Vecmem for array algebra")
add_subdirectory(algebra-plugins)

if (DETRAY_UNIT_TESTS)
  if( NOT TARGET "gtest")
    add_subdirectory(googletest)
  endif()
endif()

if (DETRAY_IO_CSV)  
  if( NOT TARGET "dfelibs" )
    add_subdirectory(dfelibs)    
  endif() 
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()
