include_directories(${MATLAB_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/toolboxes/mri_core ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math ${CMAKE_SOURCE_DIR}/toolboxes/matlab)

if (UNIX)
    if (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        set(MATLAB_SUFFIX ".mexmaci64")
    else()
        set(MATLAB_SUFFIX ".mexglnxa64")
    endif()
else()
    set(MATLAB_SUFFIX ".dll")
endif()

if (UNIX)
    add_library(gadgetron_matlab SHARED gadgetron_matlab_export.h MatlabGadget.h MatlabGadget.cpp MatlabBufferGadget.h MatlabBufferGadget.cpp )
else()
    add_library(gadgetron_matlab SHARED gadgetron_matlab_export.h MatlabGadget.h MatlabGadget.cpp)
endif()
set_target_properties(gadgetron_matlab PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
    gadgetron_matlab
    gadgetron_toolbox_log
    gadgetron_toolbox_matlab
    gadgetron_gadgetbase
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    ${MATLAB_LIBRARIES}
    ${ISMRMRD_LIBRARIES}
    optimized ${ACE_LIBRARIES}
    debug ${ACE_DEBUG_LIBRARY}
)
install(TARGETS gadgetron_matlab DESTINATION lib COMPONENT main)
install(FILES MatlabGadget.h gadgetron_matlab_export.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
install(FILES BaseGadget.m BaseBufferGadget.m bufferRecon.m scale.m accumulate_and_recon.m mask_image.m recon.m trajectoryScale.m IDEAL.m DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
install(FILES matlab.xml matlabbuffer.xml matlabnoncartesian.xml matlab_ideal_cg.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
