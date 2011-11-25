
add_definitions(-DNPY_NO_DEPRECATED_API)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${PYTHON_INCLUDE_PATH}
				       ${PYTHON_NUMPY_INCLUDE_DIR}
				       ${Boost_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_library(GadgetronPythonMRI MODULE GadgetronPythonMRI.cpp GadgetReference.cpp)
set_target_properties(GadgetronPythonMRI PROPERTIES PREFIX "")

add_library(gadgetronpython SHARED       		
	PythonGadget.cpp
	GadgetReference.cpp
	GadgetronPythonMRI.cpp
	PythonCommunicator.cpp)

target_link_libraries(gadgetronpython tinyxml optimized ${ACE_LIBRARIES}  debug ${ACE_DEBUG_LIBRARY} ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY_RELEASE})
#target_link_libraries(gadgetronpython tinyxmldebug  debug ${ACE_DEBUG_LIBRARY} ${PYTHON_DEBUG_LIBRARY} debug ${Boost_PYTHON_LIBRARY_DEBUG})

target_link_libraries(GadgetronPythonMRI tinyxml optimized ${ACE_LIBRARIES}  debug ${ACE_DEBUG_LIBRARY} optimized ${PYTHON_LIBRARY} optimized ${Boost_PYTHON_LIBRARY_RELEASE})
#target_link_libraries(GadgetronPythonMRI tinyxml debug ${ACE_DEBUG_LIBRARY} debug ${PYTHON_DEBUG_LIBRARY} debug ${Boost_PYTHON_LIBRARY_DEBUG})

install(TARGETS gadgetronpython GadgetronPythonMRI DESTINATION lib)

install(FILES 
	rms_coil_combine.py
	kspaceandimage.py
	remove_2x_oversampling.py
	accumulate_and_recon.py
	GadgetronXML.py
	image_viewer.py 
DESTINATION lib)

install(FILES
	python.xml
DESTINATION config)
