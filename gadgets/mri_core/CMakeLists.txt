IF (WIN32)
	ADD_DEFINITIONS(-DTIXML_USE_STL)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
	ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

include_directories(
  
#  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
)

find_package(Ismrmrd REQUIRED)

include_directories(  
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  )

add_library(gadgetronmricore SHARED 
  MRINoiseAdjustGadget.cpp
  PCACoilGadget.cpp
  matrix_decomposition.cpp
  matrix_vector_op.cpp
  )

target_link_libraries(gadgetronmricore tinyxml cpucore cpucore_math 
  ${ISMRMRD_LIBRARIES}  
  ${ISMRMRD_XSD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
#  ${ARMADILLO_LIBRARIES}
  )

install (FILES MRINoiseAdjustGadget.h PCACoilGadget.h DESTINATION include)

install(TARGETS gadgetronmricore DESTINATION lib)
install(FILES default_optimized.xml radial_single_optimized.xml radial_single_fixed_optimized.xml radial_ms_optimized.xml DESTINATION config)
