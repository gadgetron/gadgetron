include_directories(    ${ACE_INCLUDE_DIR} 
				${Boost_INCLUDE_DIR}  
				${CMAKE_SOURCE_DIR}/gadgets/core 
				${CMAKE_SOURCE_DIR}/apps/gadgetron
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				${CMAKE_SOURCE_DIR}/toolboxes/hostutils
				${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				${CMAKE_SOURCE_DIR}/toolboxes/hdf5utils
				${CMAKE_SOURCE_DIR}/toolboxes/hdf5mri)

add_executable(mriclient main.cpp)

target_link_libraries(mriclient tinyxml gadgettools hdf5utils hdf5mri optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS mriclient DESTINATION bin)
install (FILES ImageWriter.h HDF5ImageWriter.h DESTINATION include)
