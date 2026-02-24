# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(detray_python_bindings src/core_module_entry.cpp)

target_sources(detray_python_bindings PRIVATE src/detector.cpp)

install(TARGETS detray_python_bindings DESTINATION ${DETRAY_PYTHON_INSTALL_DIR})

set_target_properties(
    detray_python_bindings
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    detray_python_bindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/detray
)

target_link_libraries(detray_python_bindings PUBLIC detray::core_array)
