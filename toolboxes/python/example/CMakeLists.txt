include_directories(${CMAKE_SOURCE_DIR}/toolboxes/python ${ISMRMRD_INCLUDE_DIR})

add_executable(gadgetron_test_python test_python.cpp)
target_link_libraries(gadgetron_test_python gadgetron_toolbox_python)

add_executable(python_demo demo.cpp)
target_link_libraries(python_demo gadgetron_toolbox_python ${ISMRMRD_LIBRARIES})

install(TARGETS gadgetron_test_python DESTINATION bin COMPONENT main)
