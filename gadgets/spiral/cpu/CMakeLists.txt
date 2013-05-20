IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUSPIRAL__)
ENDIF (WIN32)

add_library(gadgetron_cpuspiral SHARED 
  vds.cpp 
  )

target_link_libraries(gadgetron_cpuspiral
  )

install (TARGETS gadgetron_cpuspiral DESTINATION lib)
install (FILES vds.h DESTINATION include)
