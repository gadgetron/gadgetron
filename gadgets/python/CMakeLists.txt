IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_PYTHON__)
ENDIF (WIN32)

# TODO: The Gadgetron Python code uses Numpy C-API code
# from Numpy versions < 1.7. If Numpy version is >= 1.7
# you will get compiler warnings
#add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

find_package(Ismrmrd REQUIRED)

message(STATUS ${Boost_INCLUDE_DIR} ${Boost_LIBRARIES})

include_directories(
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_BINARY_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${PYTHON_INCLUDE_PATH}
  ${NUMPY_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  )

add_library(GadgetronPythonMRI MODULE GadgetronPythonMRI.cpp GadgetReference.cpp)

#We should probably not set soversion on the python module. Causes problems with clang
#set_target_properties(GadgetronPythonMRI PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

add_library(gadgetron_python SHARED       		
	PythonGadget.cpp
	GadgetReference.cpp
	GadgetronPythonMRI.cpp
	PythonCommunicator.cpp)

set_target_properties(gadgetron_python PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_python
        gadgetron_gadgetbase
	${ISMRMRD_LIBRARIES} 
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
    ${MKL_LIBRARIES})
	
target_link_libraries(GadgetronPythonMRI
    ${ISMRMRD_LIBRARIES}
	optimized ${ACE_LIBRARIES}  
	debug ${ACE_DEBUG_LIBRARY} 
	${PYTHON_LIBRARIES} 
	${Boost_LIBRARIES}
    ${MKL_LIBRARIES})

IF (WIN32)
    SET_TARGET_PROPERTIES(GadgetronPythonMRI PROPERTIES SUFFIX .pyd)
    SET_TARGET_PROPERTIES(gadgetron_python PROPERTIES LINK_FLAGS "/LIBPATH:${PYTHON_INCLUDE_DIR}/../libs" )
ENDIF(WIN32)

set_target_properties(GadgetronPythonMRI PROPERTIES PREFIX "")

install(TARGETS gadgetron_python DESTINATION lib COMPONENT main)
install(TARGETS GadgetronPythonMRI DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(FILES 
	rms_coil_combine.py
	kspaceandimage.py
	remove_2x_oversampling.py
	accumulate_and_recon.py
	GadgetronXML.py
	image_viewer.py 
    DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(FILES
	python.xml python_short.xml
    DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
