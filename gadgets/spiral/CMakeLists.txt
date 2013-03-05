find_package(Ismrmrd REQUIRED)

include_directories(   ${ACE_INCLUDE_DIR}
                       ${CUDA_INCLUDE_DIRS}
                       ${Boost_INCLUDE_DIR}
                       ${ISMRMRD_INCLUDE_DIR}
                       ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
                       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
                       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
                       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
                       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
                       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
                       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
                       ${CMAKE_SOURCE_DIR}/toolboxes/solvers
                       ${CMAKE_SOURCE_DIR}/toolboxes/operators
                       ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
                       ${CMAKE_SOURCE_DIR}/gadgets/cgsense
		       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetronspiral SHARED SpiralGadget.cpp SpiralGadgetSW.cpp FlowPhaseSubtractionGadget.cpp vds.cpp)

target_link_libraries(gadgetronspiral gpunfft tinyxml solvers operators gpuparallelmri ${ISMRMRD_LIBRARIES} ${ISMRMRD_XSD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (TARGETS gadgetronspiral DESTINATION lib)
install (FILES spiral.xml  spiral_MSH.xml DESTINATION config)

