include_directories(${VTKWRAPPER_DIR}/ttkBottleneckDistance/)

if(UNIX AND NOT APPLE)
  ttk_add_paraview_plugin(ttkTrackingFromPersistenceDiagrams
    SOURCES ${VTKWRAPPER_DIR}/ttkTrackingFromPersistenceDiagrams/ttkTrackingFromPersistenceDiagrams.cpp
    PLUGIN_VERSION 1.0
    PLUGIN_XML TrackingFromPersistenceDiagrams.xml
    LINK trackingFromPersistenceDiagrams)
else()
  ttk_add_paraview_plugin(ttkTrackingFromPersistenceDiagrams
    SOURCES ${VTKWRAPPER_DIR}/ttkTrackingFromPersistenceDiagrams/ttkTrackingFromPersistenceDiagrams.cpp
    PLUGIN_VERSION 1.0
    PLUGIN_XML TrackingFromPersistenceDiagrams.xml
    LINK trackingFromPersistenceDiagrams ttkBottleneckDistance)
endif()
