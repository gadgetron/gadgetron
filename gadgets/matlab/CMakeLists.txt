find_package(Ismrmrd REQUIRED)
find_package(GMatlab REQUIRED)

include_directories(   ${MATLAB_INCLUDE_DIR}
					   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR} 
				       ${ISMRMRD_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

if (UNIX)
	if (APPLE)
		SET(MATLAB_SUFFIX ".mexmaci64")
	else(APPLE)
		SET(MATLAB_SUFFIX ".mexglnxa64")
	endif(APPLE)
else(UNIX)
	SET(MATLAB_SUFFIX ".dll")
endif(UNIX)

add_library(GadgetronMessage MODULE GadgetronMessage.cpp)
target_link_libraries(GadgetronMessage ${MATLAB_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})
SET_TARGET_PROPERTIES(GadgetronMessage PROPERTIES SUFFIX ${MATLAB_SUFFIX} PREFIX "")

add_library(gadgetronmatlab SHARED MatlabCommunicator.cpp MatlabGadget.cpp)
target_link_libraries(gadgetronmatlab ${MATLAB_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})
	
install(TARGETS gadgetronmatlab DESTINATION lib)
install(TARGETS GadgetronMessage DESTINATION matlab)

#install(FILES DESTINATION matlab)
install(FILES matlab.xml DESTINATION config)
