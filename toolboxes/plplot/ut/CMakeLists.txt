enable_testing()

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories( ${GTEST_INCLUDE_DIRS} 
                     ${CMAKE_SOURCE_DIR}/toolboxes/plplot 
                     ${CMAKE_SOURCE_DIR}/toolboxes/image_io )

include_directories(${GADGETRON_INCLUDE_DIR})

link_libraries( ${GTEST_LIBRARIES}
                ${ISMRMRD_LIBRARIES}
                ${ARMADILLO_LIBRARIES} 
                gadgetron_toolbox_plplot 
                gadgetron_toolbox_image_analyze_io )

add_executable(gadgetron_ut_plplot
    plplot_ut.cpp 
    plplot_test.cpp )

install(TARGETS gadgetron_ut_plplot DESTINATION bin COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)