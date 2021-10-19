add_detray_test(grids_axis
                grids_axis.cpp LINK_LIBRARIES detray::core)

add_detray_test(grids_grid2
                grids_grid2.cpp LINK_LIBRARIES detray::core)       

add_detray_test(grids_populator
                grids_populator.cpp LINK_LIBRARIES detray::core)

add_detray_test(grids_serializer
                grids_serializer.cpp LINK_LIBRARIES detray::core)
         
add_detray_test(utils_enumerate
                utils_enumerate.cpp LINK_LIBRARIES detray::core)

add_detray_test(utils_local_object_finder
                utils_local_object_finder.cpp LINK_LIBRARIES detray::core)

add_detray_test(utils_quadratic_equation
                utils_quadratic_equation.cpp LINK_LIBRARIES detray::core)
