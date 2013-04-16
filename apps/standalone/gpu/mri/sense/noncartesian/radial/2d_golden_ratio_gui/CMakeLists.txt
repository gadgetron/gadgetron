find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include( ${QT_USE_FILE} )

#We need binary and source dirs in this because of the header files created by the make process
include_directories(	${CMAKE_CURRENT_BINARY_DIR}
						${CMAKE_CURRENT_SOURCE_DIR}
						${OPENGL_INCLUDE_DIR}
						${GLUT_INCLUDE_DIR}
						${GLEW_INCLUDE_DIR}
						${CUDA_INCLUDE_DIRS}
						${Boost_INCLUDE_DIR} 
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
						${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

set(UI_UIFILES reconBaseWidget.ui radialSenseAppBaseMainWidget.ui)
qt4_wrap_ui( UI_HEADERS ${UI_UIFILES} )
 
set(UI_MOC_HEADERS radialSenseAppMainWidget.h reconWidget.h GLReconWidget.h)
qt4_wrap_cpp (UI_MOC_OUTFILES ${UI_MOC_HEADERS})

add_executable(radial_sense_gr_gui main.cpp ${UI_MOC_OUTFILES}
radialSenseAppMainWidget.cpp reconWidget.cpp GLReconWidget.cpp ${UI_HEADERS} )

target_link_libraries(radial_sense_gr_gui gpucore gpuparallelmri
gpunfft hostutils solvers ${CUDA_LIBRARIES} ${QT_QTGUI_LIBRARY} ${GLEW_LIBRARY}
${QT_QTCORE_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_gl_LIBRARY} )

if (WIN32)
set_target_properties( radial_sense_gr_gui PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE") 
endif (WIN32)

install(TARGETS radial_sense_gr_gui DESTINATION bin)
