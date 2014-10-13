find_package(Ismrmrd REQUIRED)

set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  find_package(Boost COMPONENTS thread system date_time chrono REQUIRED)
else(WIN32)
  find_package(Boost COMPONENTS thread system REQUIRED)
endif(WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
    ${CMAKE_BINARY_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools/ismrmrd
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/application
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    )

add_executable(mriclient main.cpp ImageWriter.h HDF5ImageWriter.h ImageAttribWriter.h HDF5ImageAttribWriter.h BlobFileWriter.h BlobFileWithAttribWriter.h )

add_executable(gt_alive gt_alive.cpp)

add_executable(gtdependencyquery gt_query.cpp DependencyQueryReader.h gtquery.xml)

target_link_libraries(mriclient gadgetron_toolbox_cpucore gadgetron_toolbox_gtplus_io)

target_link_libraries(gt_alive gadgetron_toolbox_cpucore 
                               gadgetron_toolbox_gadgettools 
                               optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                               ${Boost_LIBRARIES} 
                               ${ISMRMRD_LIBRARIES} )

target_link_libraries(gtdependencyquery gadgetron_toolbox_cpucore 
                                        gadgetron_toolbox_gadgettools 
                                        optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                                        ${Boost_LIBRARIES} 
                                        ${ISMRMRD_LIBRARIES} )

target_link_libraries(mriclient gadgetron_toolbox_gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(mriclient gadgetron_mricore)

install(TARGETS mriclient gt_alive gtdependencyquery DESTINATION bin COMPONENT main)
install(FILES ImageWriter.h HDF5ImageWriter.h ImageAttribWriter.h HDF5ImageAttribWriter.h BlobFileWriter.h BlobFileWithAttribWriter.h DependencyQueryReader.h DESTINATION include COMPONENT main)
install(FILES isalive.xml gtquery.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
