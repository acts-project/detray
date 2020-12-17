enable_testing()

add_detray_test(smatrix_planar_intersection
                smatrix_planar_intersection.cpp
                detray::smatrix)
                
target_link_libraries(smatrix_planar_intersection ROOT::MathCore)
