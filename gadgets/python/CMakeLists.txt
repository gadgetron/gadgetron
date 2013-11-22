IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_PYTHON__)
ENDIF (WIN32)

# All NumPy C-API code is guaranteed to compile with NumPy version >= 1.7
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${PYTHON_INCLUDE_PATH}
  ${NUMPY_INCLUDE_DIRS}
  )

add_library(GadgetronPythonMRI MODULE GadgetronPythonMRI.cpp GadgetReference.cpp)

add_library(gadgetron_python SHARED       		
	PythonGadget.cpp
	GadgetReference.cpp
	GadgetronPythonMRI.cpp
	PythonCommunicator.cpp)

target_link_libraries(gadgetron_python
	${ISMRMRD_LIBRARIES} 
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES})
	
target_link_libraries(GadgetronPythonMRI
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES} 
	${Boost_LIBRARIES})

IF (WIN32)
    SET_TARGET_PROPERTIES(GadgetronPythonMRI PROPERTIES SUFFIX .pyd)
	SET_TARGET_PROPERTIES(gadgetron_python PROPERTIES LINK_FLAGS "/LIBPATH:${PYTHON_INCLUDE_DIR}/../libs" )
ENDIF(WIN32)

set_target_properties(GadgetronPythonMRI PROPERTIES PREFIX "")
	
install(TARGETS gadgetron_python GadgetronPythonMRI DESTINATION lib)

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
