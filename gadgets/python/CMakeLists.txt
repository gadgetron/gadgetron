
add_definitions(-DNPY_NO_DEPRECATED_API)

include_directories(${ACE_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
				       ${PYTHON_INCLUDE_PATH}
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

target_link_libraries(gadgetronpython ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${PYTHON_LIBRARIES} boost_python)

target_link_libraries(GadgetronPythonMRI ${ACE_LIBRARIES} ${PYTHON_LIBRARIES}  boost_python)

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
