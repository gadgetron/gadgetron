    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_EPI__)
        link_directories(${Boost_LIBRARY_DIRS})
    endif (WIN32)

    include_directories(
        ${ACE_INCLUDE_DIR} 
        ${Boost_INCLUDE_DIR}
        ${XSD_INCLUDE_DIR}
        ${FFTW3_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${MKL_INCLUDE_DIR}
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
        ${HDF5_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIR}/cpp
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
        ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
        ${CMAKE_SOURCE_DIR}/apps/gadgetron
        ${CMAKE_SOURCE_DIR}/apps/matlab
        ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
        ${CMAKE_SOURCE_DIR}/gadgets/gtPlus 
    )

    add_library(epi ${LIBTYPE}
     EPIExport.h
     EPIReconXObject.h
     EPIReconXObjectFlat.h
     EPIReconXObjectTrapezoid.h)

    set_target_properties(epi PROPERTIES LINKER_LANGUAGE CXX)

    target_link_libraries(epi cpucore cpucore_math ${ARMADILLO_LIBRARIES})

    install (FILES 
        EPIExport.h 
        EPIReconXObject.h
        EPIReconXObjectFlat.h
        EPIReconXObjectTrapezoid.h
        DESTINATION include)

    # install(TARGETS epi DESTINATION lib)

    #if (GTEST_FOUND)
    #    add_subdirectory(ut)
    #endif (GTEST_FOUND)

