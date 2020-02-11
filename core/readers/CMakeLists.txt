
add_library(gadgetron_core_readers SHARED
        AcquisitionReader.h
        AcquisitionReader.cpp
        WaveformReader.h
        WaveformReader.cpp
        BufferReader.cpp BufferReader.h
        IsmrmrdImageArrayReader.cpp IsmrmrdImageArrayReader.h
        ImageReader.cpp ImageReader.h
        AcquisitionBucketReader.cpp AcquisitionBucketReader.h)

target_link_libraries(gadgetron_core_readers
        gadgetron_core
        ISMRMRD
        gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_core_readers PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_core_readers
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )

