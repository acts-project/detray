# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/detray)

configure_file(setup.sh.in ${_python_dir}/setup.sh @ONLY)
install(FILES ${_python_dir}/setup.sh DESTINATION python)

# Set up the init file for the detray python module
set(_target ${_python_dir}/detray/__init__.py)
file(
    CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/detray/__init__.py ${_target}
    SYMBOLIC
)

# Allow access to the python detectors package from the main __init__.py
file(
    CREATE_LINK
        ${CMAKE_CURRENT_SOURCE_DIR}/detray/detectors
        ${_python_dir}/detray/detectors
    SYMBOLIC
)

install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/detray/python/__init__.py
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}
)

# Create the metadata headers
message(STATUS "Configuring detray detector types...")
#include(FindPythonInterp)

# Generate the inbuilt metadata headers
#for each
#message(STATUS "Created xxx metadata")
#execute_process(
#   COMMAND ${PYTHON3} ${CMAKE_CURRENT_SOURCE_DIR}/my_python_script.py GET_LIST_OF_FILES
#   RESULT_VARIABLE FILES
#)
#add_custom_command(
#   OUTPUT ${FILES}
#   COMMAND ${PYTHON3} ${CMAKE_CURRENT_SOURCE_DIR}/my_python_script.py
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/my_python_script.py
#)

# Generate the custom metadata headers

# Build the python bindings
if(DETRAY_BUILD_PYTHON_BINDINGS)
    add_subdirectory(detray)
endif()
