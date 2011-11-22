IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)

include_directories(   ${ACE_INCLUDE_DIR} 
					   ${Boost_INCLUDE_DIR}  
				       ${CMAKE_SOURCE_DIR}/gadgets/core 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron
				       ${CMAKE_SOURCE_DIR}/apps/gadgetdatasender  
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_executable(mriclient main.cpp)

target_link_libraries(mriclient tinyxml gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS mriclient DESTINATION bin)
