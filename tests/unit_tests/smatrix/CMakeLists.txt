enable_testing()

add_detray_test(smatrix_masks_container
                smatrix_masks_container.cpp 
                detray::smatrix)

add_detray_test(smatrix_annulus2
                smatrix_annulus2.cpp 
                detray::smatrix)

add_detray_test(smatrix_cylinder3
                smatrix_cylinder3.cpp 
                detray::smatrix)

add_detray_test(smatrix_rectangle2
                smatrix_rectangle2.cpp 
                detray::smatrix)

add_detray_test(smatrix_ring2
                smatrix_ring2.cpp 
                detray::smatrix)

add_detray_test(smatrix_single3
                smatrix_single3.cpp 
                detray::smatrix)

add_detray_test(smatrix_trapezoid2
                smatrix_trapezoid2.cpp 
                detray::smatrix)

add_detray_test(smatrix_unmasked
                smatrix_unmasked.cpp 
                detray::smatrix)

add_detray_test(smatrix_plugin
                smatrix_plugin.cpp
                detray::smatrix)

add_detray_test(smatrix_core
                smatrix_core.cpp
                detray::smatrix)

add_detray_test(smatrix_cylinder_intersection
                smatrix_cylinder_intersection.cpp
                detray::smatrix)

add_detray_test(smatrix_planar_intersection
                smatrix_planar_intersection.cpp
                detray::smatrix)

add_detray_test(smatrix_read_detector
                smatrix_read_detector.cpp
                detray::smatrix)
                
target_link_libraries(smatrix_masks_container PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_annulus2 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_cylinder3 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_rectangle2 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_ring2 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_single3 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_trapezoid2 PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_unmasked PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_plugin PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_core PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_cylinder_intersection PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_planar_intersection PUBLIC ROOT::MathCore)
target_link_libraries(smatrix_read_detector PUBLIC ROOT::MathCore)

