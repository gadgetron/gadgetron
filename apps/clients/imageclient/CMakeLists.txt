include( ${QT_USE_FILE} ) 

#We need binary and source dirs in this because of the header files created by the make process
include_directories(    ${ACE_INCLUDE_DIR} 
			${Boost_INCLUDE_DIR}
                        ${CMAKE_CURRENT_BINARY_DIR} 
                        ${CMAKE_CURRENT_SOURCE_DIR} 
			${CMAKE_SOURCE_DIR}/gadgets/core 
			${CMAKE_SOURCE_DIR}/apps/gadgetron
			${CMAKE_SOURCE_DIR}/toolboxes/ndarray
			${CMAKE_SOURCE_DIR}/toolboxes/hostutils
			${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

set(UI_UIFILES imageClient.ui)
qt4_wrap_ui( UI_HEADERS ${UI_UIFILES} )

set(UI_MOC_HEADERS mainWindow.h)
qt4_wrap_cpp (UI_MOC_OUTFILES ${UI_MOC_HEADERS})

add_executable(imageclient main.cpp mainWindow.cpp ${UI_MOC_OUTFILES} ${UI_HEADERS})

target_link_libraries(imageclient tinyxml gadgettools 
optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTCORE_LIBRARY})

install(TARGETS imageclient DESTINATION bin)
#install (FILES ImageWriter.h DESTINATION include)
