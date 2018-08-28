if (WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif ()

include_directories( 
                    ${CMAKE_SOURCE_DIR}/toolboxes/core 
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu 
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
                    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
                    ${CMAKE_SOURCE_DIR}/toolboxes/operators
                    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
                    ${ACE_INCLUDE_DIR}
                    ${ISMRMRD_INCLUDE_DIR} )

add_executable(cpu_denoise_TV denoise_TV.cpp)
target_link_libraries(cpu_denoise_TV 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math 
                    gadgetron_toolbox_hostutils
                    ${ARMADILLO_LIBRARIES})

if (MKL_FOUND)
    target_link_libraries(cpu_denoise_TV 
                    ${MKL_LIBRARIES})
else ()
    target_link_libraries(cpu_denoise_TV 
        ${BLAS_LIBRARIES}
        ${CBLAS_LIBRARY}
    )
endif ()



install(TARGETS cpu_denoise_TV DESTINATION bin COMPONENT main)
