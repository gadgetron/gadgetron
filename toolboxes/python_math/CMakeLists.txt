find_library (BOOST_NUMPY_LIBRARY boost_numpy )

if (NOT BOOST_NUMPY_LIBRARY)
  message(FATAL_ERROR "Unable to locate boost_numpy")
endif()

message("BOOST_NUMPY_LIBRARY: ${BOOST_NUMPY_LIBRARY}")

include_directories(${ACE_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIR} 
		    ${PYTHON_INCLUDE_PATH}
		    ${NUMPY_INCLUDE_DIRS}
		    ${CMAKE_SOURCE_DIR}/toolboxes/core
		    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                    )

add_library(gadgetron_toolbox_python_math SHARED
  PythonMath.cpp
  PythonMath.h
  python_math_export.h
)

target_link_libraries(gadgetron_toolbox_python_math
                      gadgetron_toolbox_cpucore
                      	${PYTHON_LIBRARIES}
			${Boost_LIBRARIES}
			${BOOST_NUMPY_LIBRARY}
			${ACE_LIBRARIES}
			gadgetron_toolbox_log)

add_executable(gadgetron_test_python_math
  test_python_math.cpp)

target_link_libraries(gadgetron_test_python_math
  gadgetron_toolbox_python_math
  gadgetron_toolbox_cpucore
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES}
  ${BOOST_NUMPY_LIBRARY}
  gadgetron_toolbox_log)

set_target_properties(gadgetron_toolbox_python_math  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_toolbox_python_math DESTINATION lib COMPONENT main)
install(TARGETS gadgetron_test_python_math DESTINATION bin COMPONENT main)

install(FILES 
  PythonMath.h
  python_math_export.h
  DESTINATION include COMPONENT main)
