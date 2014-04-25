IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CARTESIAN__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${ISMRMRD_XSD_INCLUDE_DIR}
  )

add_library(gadgetron_cartesian SHARED 
  gadgetron_cartesian_export.h
  CartesianToGenericGadget.h
  CartesianToGenericGadget.cpp
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_cartesian cpucore 
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_cartesian cpucore_math )
endif(ARMADILLO_FOUND)

install (FILES  gadgetron_cartesian_export.h
                CartesianToGenericGadget.h
                DESTINATION include)

install (TARGETS gadgetron_cartesian DESTINATION lib)
