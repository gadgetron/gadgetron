include_directories(${ISMRMRD_INCLUDE_DIR})

add_executable(mrd_to_stream mrd_to_stream.cpp)
add_executable(stream_to_mrd stream_to_mrd.cpp)


target_include_directories(mrd_to_stream
    PRIVATE
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(stream_to_mrd
    PRIVATE
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_BINARY_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image

    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})


target_link_libraries(mrd_to_stream
    Boost::program_options
    Boost::system
    Boost::filesystem
    ISMRMRD::ISMRMRD)

target_link_libraries(stream_to_mrd
    Boost::program_options
    Boost::system
    Boost::filesystem
    ISMRMRD::ISMRMRD
    gadgetron_core_readers)


install(TARGETS stream_to_mrd DESTINATION bin COMPONENT main)
install(TARGETS mrd_to_stream DESTINATION bin COMPONENT main)
