include_directories(${ACE_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/gadgets/core 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron
				       ${CMAKE_SOURCE_DIR}/apps/gadgetdatasender  
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_executable(mriclient main.cpp)

target_link_libraries(mriclient ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} gadgettools)

install(TARGETS mriclient DESTINATION bin)
