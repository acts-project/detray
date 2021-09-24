enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(vc_array_${etest}
      vc_array_${etest}.cpp
      LINK_LIBRARIES
      detray::vc_array)
endforeach(etest)

if (DETRAY_IO_CSV)
    add_detray_test(vc_array_read_detector
      vc_array_read_detector.cpp
      LINK_LIBRARIES
      detray::vc_array)
endif()

if (DETRAY_DISPLAY)
    add_detray_test(vc_array_display_masks
      vc_array_display_masks.cpp
      LINK_LIBRARIES
      detray::vc_array detray::display)
endif()
