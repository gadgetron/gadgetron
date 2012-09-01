find_package(Ismrmrd REQUIRED)
find_package(HDF5 1.8 COMPONENTS C CXX REQUIRED)

include_directories(    
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}  
  ${HDF5_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}/cpp
  ${CMAKE_SOURCE_DIR}/gadgets/core 
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
  ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${ISMRMRD_INCLUDE_DIR}
  )

add_executable(mriclient main.cpp)

target_link_libraries(mriclient tinyxml gadgettools  ${HDF5_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ISMRMRD_LIBRARIES})


install(TARGETS mriclient DESTINATION bin)
install (FILES ImageWriter.h HDF5ImageWriter.h DESTINATION include)
install(FILES ${ISMRMRD_LIBRARIES} DESTINATION lib)
