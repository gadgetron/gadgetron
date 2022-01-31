enable_testing()

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif()

message("Found PLPLOT, build plplot test ...")
add_executable(gadgetron_ut_plplot
    plplot_ut.cpp 
    plplot_test.cpp )

target_link_libraries(gadgetron_ut_plplot gadgetron_toolbox_plplot gadgetron_toolbox_image_analyze_io ${GTEST_LIBRARIES})
install(TARGETS gadgetron_ut_plplot DESTINATION bin COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)