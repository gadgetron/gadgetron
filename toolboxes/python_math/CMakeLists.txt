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
			${ACE_LIBRARIES}
			gadgetron_toolbox_log)

set_target_properties(gadgetron_toolbox_python_math  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties (gadgetron_toolbox_python_math PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_TOOLBOX_PYTHON_MATH__")

install(TARGETS gadgetron_toolbox_python_math DESTINATION lib COMPONENT main)

install(FILES 
  PythonMath.h
  python_math_export.h
  DESTINATION include COMPONENT main)

add_subdirectory(example)
