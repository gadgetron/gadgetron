
add_library(gadgetron_core_writers SHARED
        ImageWriter.cpp
        ImageWriter.h
        BufferWriter.cpp
        BufferWriter.h
        IsmrmrdImageArrayWriter.cpp
        IsmrmrdImageArrayWriter.h
        AcquisitionBucketWriter.cpp
        AcquisitionBucketWriter.h
        AcquisitionWriter.cpp
        AcquisitionWriter.h
        WaveformWriter.cpp
        WaveformWriter.h)

target_link_libraries(gadgetron_core_writers
        gadgetron_core
        ISMRMRD::ISMRMRD
        gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_core_writers PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_core_writers
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)