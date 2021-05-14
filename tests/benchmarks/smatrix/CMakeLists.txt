add_detray_benchmark(smatrix_intersect_all
                     smatrix_intersect_all.cpp
                     algebra::smatrix)

add_detray_benchmark(smatrix_intersect_surfaces
                     smatrix_intersect_surfaces.cpp
                     algebra::smatrix)

add_detray_benchmark(smatrix_masks
                     smatrix_masks.cpp
                     algebra::smatrix)

target_link_libraries(smatrix_intersect_all ROOT::MathCore)
target_link_libraries(smatrix_intersect_surfaces ROOT::MathCore)
target_link_libraries(smatrix_masks ROOT::MathCore)
