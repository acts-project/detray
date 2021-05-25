message(STATUS "Building 'detray::smatrix' plugin")

add_library(detray_smatrix INTERFACE)

target_include_directories(detray_smatrix
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_smatrix
    INTERFACE -DALGEBRA_CUSTOM_SCALARTYPE=${ALGEBRA_CUSTOM_SCALARTYPE})
endif()

target_link_libraries(detray_smatrix INTERFACE algebra::smatrix)

add_library(detray::smatrix ALIAS detray_smatrix)
