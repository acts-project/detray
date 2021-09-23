add_subdirectory(algebra-plugins)

if (DETRAY_UNIT_TESTS)
  if (NOT gtest)
    add_subdirectory(googletest)
  endif()
endif()

if (DETRAY_IO_CSV)
  if( TARGET "dfelibs" )
    return()
  elseif()
    add_subdirectory(dfelibs)    
  endif() 
endif()

if(DETRAY_DISPLAY)
    add_subdirectory(matplotplusplus)
endif()
