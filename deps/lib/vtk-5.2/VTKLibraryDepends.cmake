# Generated by CMake 2.6-patch 0

IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
  # Information for CMake 2.6 and above.
  SET("QVTKWidgetPlugin_LIB_DEPENDS" "optimized;C:/Qt/4.4/lib/QtCore4.lib;debug;C:/Qt/4.4/lib/QtCored4.lib;optimized;C:/Qt/4.4/lib/QtGui4.lib;debug;C:/Qt/4.4/lib/QtGuid4.lib;")
  SET("QVTK_LIB_DEPENDS" "optimized;C:/Qt/4.4/lib/QtGui4.lib;debug;C:/Qt/4.4/lib/QtGuid4.lib;general;imm32;general;winmm;optimized;C:/Qt/4.4/lib/QtSql4.lib;debug;C:/Qt/4.4/lib/QtSqld4.lib;optimized;C:/Qt/4.4/lib/QtCore4.lib;debug;C:/Qt/4.4/lib/QtCored4.lib;general;ws2_32;general;vtkRendering;general;vtkGraphics;general;vtkImaging;general;vtkCommon;general;vtkViews;")
  SET("vtkCommon_LIB_DEPENDS" "general;vtksys;general;wsock32;")
  SET("vtkFiltering_LIB_DEPENDS" "general;vtkCommon;")
  SET("vtkGenericFiltering_LIB_DEPENDS" "general;vtkFiltering;general;vtkGraphics;")
  SET("vtkGraphics_LIB_DEPENDS" "general;vtkFiltering;general;vtkverdict;")
  SET("vtkHybrid_LIB_DEPENDS" "general;vtkRendering;general;vtkIO;general;vtkexoIIc;general;vfw32;")
  SET("vtkIO_LIB_DEPENDS" "general;vtkFiltering;general;vtkDICOMParser;general;vtkNetCDF;general;vtkmetaio;general;vtksqlite;general;vtkpng;general;vtkzlib;general;vtkjpeg;general;vtktiff;general;vtkexpat;general;vtksys;general;vfw32;")
  SET("vtkImaging_LIB_DEPENDS" "general;vtkFiltering;")
  SET("vtkInfovis_LIB_DEPENDS" "general;vtkWidgets;general;vtklibxml2;")
  SET("vtkRendering_LIB_DEPENDS" "general;vtkGraphics;general;vtkImaging;general;vtkIO;general;vtkftgl;general;vtkfreetype;general;vtkzlib;general;vtkpng;general;opengl32;")
  SET("vtkViews_LIB_DEPENDS" "general;vtkInfovis;")
  SET("vtkVolumeRendering_LIB_DEPENDS" "general;vtkRendering;general;vtkIO;")
  SET("vtkWidgets_LIB_DEPENDS" "general;vtkRendering;general;vtkHybrid;")
  SET("vtkexoIIc_LIB_DEPENDS" "general;vtkNetCDF;")
  SET("vtkftgl_LIB_DEPENDS" "general;opengl32;general;vtkfreetype;")
  SET("vtklibxml2_LIB_DEPENDS" "general;vtkzlib;")
  SET("vtkmetaio_LIB_DEPENDS" "general;vtkzlib;general;vtksys;general;comctl32;general;wsock32;")
  SET("vtkpng_LIB_DEPENDS" "general;vtkzlib;")
  SET("vtksys_LIB_DEPENDS" "general;ws2_32;")
  SET("vtktiff_LIB_DEPENDS" "general;vtkzlib;general;vtkjpeg;")
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
  # Information for CMake 2.4 and lower.
  SET("QVTKWidgetPlugin_LIB_DEPENDS" "C:/Qt/4.4/lib/QtCore4.lib;C:/Qt/4.4/lib/QtCored4.lib;C:/Qt/4.4/lib/QtGui4.lib;C:/Qt/4.4/lib/QtGuid4.lib;")
  SET("QVTK_LIB_DEPENDS" "C:/Qt/4.4/lib/QtGui4.lib;C:/Qt/4.4/lib/QtGuid4.lib;imm32;winmm;C:/Qt/4.4/lib/QtSql4.lib;C:/Qt/4.4/lib/QtSqld4.lib;C:/Qt/4.4/lib/QtCore4.lib;C:/Qt/4.4/lib/QtCored4.lib;ws2_32;vtkRendering;vtkGraphics;vtkImaging;vtkCommon;vtkViews;")
  SET("vtkCommon_LIB_DEPENDS" "vtksys;wsock32;")
  SET("vtkFiltering_LIB_DEPENDS" "vtkCommon;")
  SET("vtkGenericFiltering_LIB_DEPENDS" "vtkFiltering;vtkGraphics;")
  SET("vtkGraphics_LIB_DEPENDS" "vtkFiltering;vtkverdict;")
  SET("vtkHybrid_LIB_DEPENDS" "vtkRendering;vtkIO;vtkexoIIc;vfw32;")
  SET("vtkIO_LIB_DEPENDS" "vtkFiltering;vtkDICOMParser;vtkNetCDF;vtkmetaio;vtksqlite;vtkpng;vtkzlib;vtkjpeg;vtktiff;vtkexpat;vtksys;vfw32;")
  SET("vtkImaging_LIB_DEPENDS" "vtkFiltering;")
  SET("vtkInfovis_LIB_DEPENDS" "vtkWidgets;vtklibxml2;")
  SET("vtkRendering_LIB_DEPENDS" "vtkGraphics;vtkImaging;vtkIO;vtkftgl;vtkfreetype;vtkzlib;vtkpng;opengl32;")
  SET("vtkViews_LIB_DEPENDS" "vtkInfovis;")
  SET("vtkVolumeRendering_LIB_DEPENDS" "vtkRendering;vtkIO;")
  SET("vtkWidgets_LIB_DEPENDS" "vtkRendering;vtkHybrid;")
  SET("vtkexoIIc_LIB_DEPENDS" "vtkNetCDF;")
  SET("vtkftgl_LIB_DEPENDS" "opengl32;vtkfreetype;")
  SET("vtklibxml2_LIB_DEPENDS" "vtkzlib;")
  SET("vtkmetaio_LIB_DEPENDS" "vtkzlib;vtksys;comctl32;wsock32;")
  SET("vtkpng_LIB_DEPENDS" "vtkzlib;")
  SET("vtksys_LIB_DEPENDS" "ws2_32;")
  SET("vtktiff_LIB_DEPENDS" "vtkzlib;vtkjpeg;")
  SET("C:/Qt/4.4/lib/QtCore4.lib_LINK_TYPE" "optimized")
  SET("C:/Qt/4.4/lib/QtCored4.lib_LINK_TYPE" "debug")
  SET("C:/Qt/4.4/lib/QtGui4.lib_LINK_TYPE" "optimized")
  SET("C:/Qt/4.4/lib/QtGuid4.lib_LINK_TYPE" "debug")
  SET("C:/Qt/4.4/lib/QtSql4.lib_LINK_TYPE" "optimized")
  SET("C:/Qt/4.4/lib/QtSqld4.lib_LINK_TYPE" "debug")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
