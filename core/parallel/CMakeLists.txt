
add_library(gadgetron_core_parallel SHARED
        Branch.h
		Branch.hpp
		Branch.cpp
		Merge.h
		Merge.cpp
		Fanout.h
		Fanout.hpp
		Fanout.cpp
		UnorderedMerge.h
		UnorderedMerge.cpp)

target_link_libraries(gadgetron_core_parallel
        gadgetron_core)

set_target_properties(gadgetron_core_parallel PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_core_parallel
		EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        Branch.h
        Branch.hpp
        Merge.h
        Fanout.h
        Fanout.hpp
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
