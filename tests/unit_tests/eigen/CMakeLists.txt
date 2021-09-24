enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(eigen_${etest}
                eigen_${etest}.cpp LINK_LIBRARIES detray::eigen)
endforeach(etest)

if (DETRAY_IO_CSV)
    add_detray_test(eigen_read_detector
      eigen_read_detector.cpp
      LINK_LIBRARIES
      detray::eigen)
endif()

if (DETRAY_DISPLAY)
    add_detray_test(eigen_display_masks
      eigen_display_masks.cpp
      LINK_LIBRARIES
      detray::eigen detray::display)
endif()
