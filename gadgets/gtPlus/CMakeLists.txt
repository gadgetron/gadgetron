
include_directories( 
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/toolboxes/cloudbus
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools/ismrmrd
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
    ${CMAKE_SOURCE_DIR}/gadgets/gtPlus 
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
  )

IF (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
    ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)

set( gtCloud_files 
        config/gtCloud/myCloud_2DT.txt 
        config/gtCloud/myCloud_3DT.txt 
        config/gtCloud/myCloud_2DT_DualLayer.txt 
        config/gtCloud/myCloud_2DT_DualLayer_FirstLayer.txt )

set( config_files 

    config/GT_2DT_Cartesian.xml
    config/GT_2DT_Cartesian_CloudNode.xml
    config/GT_2DT_Cartesian_Dicom.xml
    config/GT_2DT_Cartesian_DualLayer_Gateway_L1SPIRIT.xml
    config/GT_2DT_Cartesian_DualLayer_Gateway_SPIRIT.xml
    config/GT_2DT_Cartesian_FirstLayer_CloudNode.xml
    config/GT_2DT_Cartesian_GFactor.xml
    config/GT_2DT_Cartesian_ImageTrigger_Dicom.xml
    config/GT_2DT_Cartesian_L1SPIRIT.xml
    config/GT_2DT_Cartesian_PseudoReplica_SNRUnitRecon.xml
    config/GT_2DT_Cartesian_SingleLayer_CloudNode.xml
    config/GT_2DT_Cartesian_SPIRIT.xml
    config/GT_2DT_FatWater.xml

    config/GT_2DT_HASTE.xml
    config/GT_2DT_HASTE_MOCO_AVE.xml

    config/GT_2DT_T2W.xml

    config/GT_2DT_LGE.xml

    config/GT_2DT_MOLLI.xml
    config/GT_2DT_MOLLI_Offline.xml

    config/GT_2DT_Perfusion.xml

    config/GT_2DT_PseudoReplica_SNRUnitRecon_DataExport.xml

    config/GT_2DT_RealTimeCine.xml
    config/GT_2DT_RealTimeFlow.xml

    config/GT_2DT_RTCine_L1SPIRIT_PhysioInterp.xml
    config/GT_2DT_RTCine_L1SPIRIT_PhysioInterp_DualLayer_Gateway.xml

    config/GT_3DT_Cartesian.xml
    config/GT_3DT_Cartesian_CloudNode.xml
    config/GT_3DT_Cartesian_GFactor.xml
    config/GT_3DT_Cartesian_L1SPIRIT.xml
    config/GT_3DT_Cartesian_SingleLayer_L1SPIRIT.xml )

set( gadgetronPlus_header_files GtPlusGadgetImageArray.h
                                GtPlusAccumulatorWorkOrderTriggerGadget.h
                                GtPlusAccumulatorImageTriggerGadget.h
                                GtPlusGadgetOpenMP.h
                                GtPlusReconGadget.h
                                GtPlusRecon2DTGadget.h
                                GtPlusRecon3DTGadget.h
                                GtPlusRecon2DTGadgetCloud.h
                                GtPlusRecon2DTCloudPackage.h
                                GadgetCloudJobMessageReadWrite.h
                                GtPlusReconJob2DTGadget.h 
                                GtPlusReconJob3DTGadget.h 
                                GtPlusReconJob2DTGadgetCloud.h 
                                GtPlusImageReconGadget.h 
                                GtPlusReconGadgetUtil.h 
                                )

set( gadgetronPlus_src_files GtPlusGadgetImageArray.cpp
                            GtPlusAccumulatorWorkOrderTriggerGadget.cpp
                            GtPlusAccumulatorImageTriggerGadget.cpp
                            GtPlusGadgetOpenMP.cpp
                            GtPlusReconGadget.cpp
                            GtPlusRecon2DTGadget.cpp
                            GtPlusRecon3DTGadget.cpp
                            GtPlusRecon2DTGadgetCloud.cpp
                            GadgetCloudJobMessageReadWrite.cpp
                            GtPlusReconJob2DTGadget.cpp 
                            GtPlusReconJob3DTGadget.cpp 
                            GtPlusReconJob2DTGadgetCloud.cpp 
                            GtPlusImageReconGadget.cpp 
                            GtPlusReconGadgetUtil.cpp 
                            )

set( config_gtCloud_files ${config_files} ${gtCloud_files} )
source_group(config FILES ${config_gtCloud_files})

add_library(gadgetronPlus SHARED
            GtPlusGadgetExport.h 
            ${gadgetronPlus_header_files} 
            ${gadgetronPlus_src_files} 
            ${config_gtCloud_files} 
        )

set_target_properties(gadgetronPlus PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetronPlus 
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore 
    gadgetron_toolbox_cpucore_math 
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_gtplus 
    gadgetron_toolbox_gadgettools 
    gadgetron_toolbox_cloudbus 
    gadgetron_toolbox_mri_core 
    ${Boost_LIBRARIES}
    ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    )

if (CUDA_FOUND)  
  include_directories( 
    ${CUDA_INCLUDE_DIRS}
    )
  target_link_libraries(gadgetronPlus gadgetron_toolbox_gpuparallelmri)
endif(CUDA_FOUND)

install (FILES 
        GtPlusGadgetExport.h 
        ${gadgetronPlus_header_files} 
        DESTINATION include COMPONENT main)

install (FILES  ${config_files} 
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install (FILES  ${gtCloud_files} 
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH}/gtCloud COMPONENT main)

install(TARGETS gadgetronPlus DESTINATION lib COMPONENT main)
