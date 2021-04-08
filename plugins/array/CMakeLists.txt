message(STATUS "Building 'detray::array' plugin")

add_library(detray_array INTERFACE)

target_include_directories(detray_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_array
    INTERFACE -DDETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE})
endif()

add_library(detray::array ALIAS detray_array)
