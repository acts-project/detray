message(STATUS "Building 'detray::eigen' plugin")

add_library(detray_eigen INTERFACE)

target_include_directories(detray_eigen
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_eigen
    INTERFACE -DALGEBRA_CUSTOM_SCALARTYPE=${ALGEBRA_CUSTOM_SCALARTYPE})
endif()

target_link_libraries(detray_eigen INTERFACE algebra::eigen_eigen vecmem::core)

add_library(detray::eigen ALIAS detray_eigen)
