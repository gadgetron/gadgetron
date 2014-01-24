
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
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/gadgets/core
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
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

add_library(gadgetronPlus SHARED
            GtPlusGadgetExport.h
            GtPlusGadgetImageArray.h
            GtPlusGadgetImageArray.cpp
            GtPlusAccumulatorWorkOrderTriggerGadget.h
            GtPlusAccumulatorWorkOrderTriggerGadget.cpp
            GtPlusGadgetOpenMP.h
            GtPlusGadgetOpenMP.cpp
            GtPlusReconGadget.h
            GtPlusReconGadget.cpp
            GtPlusRecon2DTGadget.h
            GtPlusRecon2DTGadget.cpp
            GtPlusRecon3DTGadget.h
            GtPlusRecon3DTGadget.cpp
            GtPlusRecon2DTGadgetCloud.h
            GtPlusRecon2DTGadgetCloud.cpp
            GtPlusRecon2DTCloudPackage.h
            GadgetCloudJobMessageReadWrite.h
            GadgetCloudJobMessageReadWrite.cpp
            GtPlusReconJob2DTGadget.h 
            GtPlusReconJob2DTGadget.cpp 
            GtPlusReconJob3DTGadget.h 
            GtPlusReconJob3DTGadget.cpp 
            GtPlusReconJob2DTGadgetCloud.h 
            GtPlusReconJob2DTGadgetCloud.cpp 
            ${ISMRMRD_XSD_SOURCE} )

target_link_libraries(gadgetronPlus 
    cpucore 
    cpucore_math 
    gtplus 
    gadgettools 
    ${MKL_LIBRARIES} 
    ${Boost_LIBRARIES}
    ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
    ${XERCESC_LIBRARIES} )

    if (CUDA_FOUND)
        target_link_libraries(gadgetronPlus gpuparallelmri)
    endif (CUDA_FOUND)

install (FILES 
        GtPlusGadgetExport.h
        DESTINATION include)

install(TARGETS gadgetronPlus DESTINATION lib)
# install(FILES default.xml default_short.xml DESTINATION config)
