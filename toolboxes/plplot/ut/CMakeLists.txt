ENABLE_TESTING()

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( ${GTEST_INCLUDE_DIRS} 
                     ${CMAKE_SOURCE_DIR}/toolboxes/plplot 
                     ${CMAKE_SOURCE_DIR}/toolboxes/gtplus 
                     ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util )

include_directories(${GADGETRON_INCLUDE_DIR})

link_libraries( ${GTEST_LIBRARIES} 
                ${ISMRMRD_LIBRARIES} 
                ${ARMADILLO_LIBRARIES} 
                gadgetron_toolbox_plplot 
                gadgetron_toolbox_gtplus_io )

add_executable(gadgetron_ut_plplot
    plplot_ut.cpp 
    plplot_test.cpp )

install(TARGETS gadgetron_ut_plplot DESTINATION bin COMPONENT main)
