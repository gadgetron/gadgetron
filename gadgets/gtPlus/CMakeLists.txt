
include_directories( 
    ${CMAKE_SOURCE_DIR}/gadgets/core
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${ISMRMRD_SCHEMA_DIR}
    ${ISMRMRD_XSD_INCLUDE_DIR}
    ${XSD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/gadgets/core
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/application
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools/ismrmrd
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
    ${CMAKE_SOURCE_DIR}/gadgets/gtPlus 
  )

IF (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
    ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

set( gadgetronPlus_header_files GtPlusGadgetImageArray.h
                                GtPlusAccumulatorWorkOrderTriggerGadget.h
                                GtPlusGadgetOpenMP.h
                                GtPlusReconGadget.h
                                GtPlusRecon2DTGadget.h
                                GtPlusRecon3DTGadget.h
                                GtPlusRecon2DTGadgetCloud.h
                                GtPlusRecon2DTCloudPackage.h
                                GadgetCloudJobMessageReadWrite.h
                                GtPlusReconJob2DTGadget.h 
                                GtPlusReconJob3DTGadget.h 
                                GtPlusReconJob2DTGadgetCloud.h )

set( gadgetronPlus_src_files GtPlusGadgetImageArray.cpp
                            GtPlusAccumulatorWorkOrderTriggerGadget.cpp
                            GtPlusGadgetOpenMP.cpp
                            GtPlusReconGadget.cpp
                            GtPlusRecon2DTGadget.cpp
                            GtPlusRecon3DTGadget.cpp
                            GtPlusRecon2DTGadgetCloud.cpp
                            GadgetCloudJobMessageReadWrite.cpp
                            GtPlusReconJob2DTGadget.cpp 
                            GtPlusReconJob3DTGadget.cpp 
                            GtPlusReconJob2DTGadgetCloud.cpp )

add_library(gadgetronPlus SHARED
            GtPlusGadgetExport.h 
            ${gadgetronPlus_header_files} 
            ${gadgetronPlus_src_files} 
            ${ISMRMRD_XSD_SOURCE} )

target_link_libraries(gadgetronPlus 
    cpucore 
    cpucore_math 
    cpufft
    gtplus 
    gadgettools 
    ${Boost_LIBRARIES}
    ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
    ${XERCESC_LIBRARIES} )
  
target_link_libraries(gadgetronPlus gpuparallelmri)


install (FILES 
        GtPlusGadgetExport.h 
        ${gadgetronPlus_header_files} 
        DESTINATION include)

install(TARGETS gadgetronPlus DESTINATION lib)
