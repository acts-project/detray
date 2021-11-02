

add_subdirectory(array)

if (DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()

if (DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()

if (DETRAY_VC_PLUGIN)
    add_subdirectory(vc)
endif()
