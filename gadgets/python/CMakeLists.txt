
add_definitions(-DNPY_NO_DEPRECATED_API)

include_directories(${ACE_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
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

target_link_libraries(gadgetronpython optimized ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY_RELEASE})
target_link_libraries(gadgetronpython debug ${ACE_DEBUG_LIBRARY} ${TINYXML_LIBRARIES} ${PYTHON_DEBUG_LIBRARY} ${Boost_PYTHON_LIBRARY_DEBUG})


target_link_libraries(GadgetronPythonMRI optimized ${ACE_LIBRARIES} ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY_RELEASE})
target_link_libraries(GadgetronPythonMRI debug ${ACE_DEBUG_LIBRARY} ${PYTHON_DEBUG_LIBRARY} ${Boost_PYTHON_LIBRARY_DEBUG})


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
