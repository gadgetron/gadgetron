if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_PR__)
endif ()


set(pr_header_fiels pr_export.h 
                    pr_kmeans.h )

set(pr_src_fiels pr_kmeans.cpp )

add_library(gadgetron_toolbox_pr SHARED 
            ${pr_header_fiels} 
            ${pr_src_fiels} )

set_target_properties(gadgetron_toolbox_pr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_pr
        			ISMRMRD::ISMRMRD
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_image_analyze_io 
                    gadgetron_toolbox_cpudwt 
                    gadgetron_toolbox_cpu_image )
target_include_directories(gadgetron_toolbox_pr
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS gadgetron_toolbox_pr
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES ${pr_header_fiels} 
      DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
