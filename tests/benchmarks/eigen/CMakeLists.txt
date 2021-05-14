add_detray_benchmark(eigen_find_volume
                     eigen_find_volume.cpp
                     detray::eigen)

add_detray_benchmark(eigen_masks
                     eigen_masks.cpp
                     algebra::eigen)

add_detray_benchmark(eigen_intersect_all
                     eigen_intersect_all.cpp
                     algebra::eigen)

add_detray_benchmark(eigen_intersect_surfaces
                     eigen_intersect_surfaces.cpp
                     algebra::eigen)
