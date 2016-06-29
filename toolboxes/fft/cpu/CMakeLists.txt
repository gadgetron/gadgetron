if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUFFT__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${FFTW3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
  )

add_library(gadgetron_toolbox_cpufft SHARED 
  cpufft_export.h 
  hoNDFFT.h
  hoNDFFT.cpp
  )

set_target_properties(gadgetron_toolbox_cpufft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

if (MKL_FOUND)
    target_link_libraries(gadgetron_toolbox_cpufft
      gadgetron_toolbox_cpucore
      gadgetron_toolbox_log
      gadgetron_toolbox_cpucore_math
      ${MKL_LIBRARIES} 
      ${Boost_LIBRARIES} 
      )
else ()
target_link_libraries(gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore_math
  ${FFTW3_LIBRARIES} 
  )
endif ()

install(TARGETS gadgetron_toolbox_cpufft DESTINATION lib COMPONENT main)

install(FILES
  cpufft_export.h 
  hoNDFFT.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
