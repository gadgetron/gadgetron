include( ${QT_USE_FILE} )

#We need binary and source dirs in this because of the header files created by the make process
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
)

set(UI_UIFILES reconBaseWidget.ui radialSenseAppBaseMainWidget.ui)
qt4_wrap_ui( UI_HEADERS ${UI_UIFILES} )
 
set(UI_MOC_HEADERS radialSenseAppMainWidget.h reconWidget.h GLReconWidget.h)
qt4_wrap_cpp (UI_MOC_OUTFILES ${UI_MOC_HEADERS})

add_executable(sense_2d_golden_radial_gui main.cpp ${UI_MOC_OUTFILES}
radialSenseAppMainWidget.cpp reconWidget.cpp GLReconWidget.cpp ${UI_HEADERS} )

target_link_libraries(sense_2d_golden_radial_gui gadgetron_toolbox_gpucore gadgetron_toolbox_gpuparallelmri
gadgetron_toolbox_gpunfft gadgetron_toolbox_hostutils gadgetron_toolbox_gpusolvers gpuoperators ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES} ${QT_QTGUI_LIBRARY} ${GLEW_LIBRARY}
${QT_QTCORE_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_gl_LIBRARY} )

if (WIN32)
set_target_properties( sense_2d_golden_radial_gui PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE") 
endif ()

install(TARGETS sense_2d_golden_radial_gui DESTINATION bin COMPONENT main)
