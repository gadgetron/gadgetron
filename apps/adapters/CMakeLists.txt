include_directories(${ISMRMRD_INCLUDE_DIR})

add_executable(hdf5_to_stream hdf5_to_stream.cpp)
add_executable(stream_to_hdf5 stream_to_hdf5.cpp)


target_include_directories(hdf5_to_stream
    PRIVATE
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(stream_to_hdf5
    PRIVATE
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_BINARY_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image

    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})


target_link_libraries(hdf5_to_stream
    Boost::program_options
    Boost::system
    Boost::filesystem
    ISMRMRD::ISMRMRD)

target_link_libraries(stream_to_hdf5
    Boost::program_options
    Boost::system
    Boost::filesystem
    ISMRMRD::ISMRMRD
    gadgetron_core_readers)


install(TARGETS stream_to_hdf5 DESTINATION bin COMPONENT main)
install(TARGETS hdf5_to_stream DESTINATION bin COMPONENT main)
