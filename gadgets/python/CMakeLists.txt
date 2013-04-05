# TODO: The Gadgetron Python code uses Numpy C-API code
# from Numpy versions < 1.7. If Numpy version is >= 1.7
# you will get compiler warnings
#add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

find_package(Ismrmrd REQUIRED)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${PYTHON_INCLUDE_PATH}
				       ${PYTHON_NUMPY_INCLUDE_DIR}
				       ${Boost_INCLUDE_DIR} 
				       ${ISMRMRD_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_library(GadgetronPythonMRI MODULE GadgetronPythonMRI.cpp GadgetReference.cpp)


add_library(gadgetronpython SHARED       		
	PythonGadget.cpp
	GadgetReference.cpp
	GadgetronPythonMRI.cpp
	PythonCommunicator.cpp)

target_link_libraries(gadgetronpython tinyxml
	${ISMRMRD_LIBRARIES} 
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES})
	
target_link_libraries(GadgetronPythonMRI tinyxml 
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES} 
	${Boost_LIBRARIES})

IF (WIN32)
    SET_TARGET_PROPERTIES(GadgetronPythonMRI PROPERTIES SUFFIX .pyd)
	SET_TARGET_PROPERTIES(gadgetronpython PROPERTIES LINK_FLAGS "/LIBPATH:${PYTHON_INCLUDE_DIR}/../libs" )
ENDIF(WIN32)

set_target_properties(GadgetronPythonMRI PROPERTIES PREFIX "")
	
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
	python.xml python_short.xml
DESTINATION config)
