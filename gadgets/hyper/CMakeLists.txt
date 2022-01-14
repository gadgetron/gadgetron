if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_HYPER__)
endif ()


add_library(gadgetron_hyper SHARED 
  NFFT2DGadget.h NFFT2DGadget.cpp
  CMRT3DGadget.h CMRT3DGadget.cpp
  CMRTGadget.h CMRTGadget.cpp
  CSIGadget.h CSIGadget.cpp
  gpuCSICoilEstimationGadget.cpp
  
  )

set_target_properties(gadgetron_hyper PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_hyper
  gadgetron_core
  gadgetron_mricore
  gadgetron_gpuparallelmri
  gadgetron_toolbox_gpu
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math)


install(FILES 
  NFFT2DGadget.h
  CMRTGadget.h
  CMRT3DGadget.h
  CSIGadget.h
  gadgetron_hyper_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH})

install(FILES NFFT2D.xml CMRT3D.xml CMRT.xml FS-CSI.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})


install(TARGETS gadgetron_hyper
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)