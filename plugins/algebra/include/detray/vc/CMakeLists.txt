message(STATUS "Building 'detray::vc_array' plugin")

add_library(detray_vc_array INTERFACE)

target_include_directories(detray_vc_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_vc_array
    INTERFACE -DALGEBRA_CUSTOM_SCALARTYPE=${ALGEBRA_CUSTOM_SCALARTYPE})
endif()

target_link_libraries(detray_vc_array INTERFACE algebra::vc_array Vc)

add_library(detray::vc_array ALIAS detray_vc_array)
