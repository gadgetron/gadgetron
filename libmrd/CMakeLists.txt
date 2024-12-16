include(FetchContent)
FetchContent_Declare(mrd2
    GIT_REPOSITORY https://github.com/ismrmrd/mrd.git
    GIT_TAG v2.1.0
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(mrd2)

set(HOWARD_HINNANT_DATE_MINIMUM_VERSION "3.0.0")
find_package(date ${HOWARD_HINNANT_DATE_MINIMUM_VERSION} REQUIRED)

set(XTENSOR_MINIMUM_VERSION "0.21.10")
find_package(xtensor ${XTENSOR_MINIMUM_VERSION} REQUIRED)

set(NLOHMANN_JSON_MINIMUM_VERSION "3.11.1")
find_package(nlohmann_json ${NLOHMANN_JSON_MINIMUM_VERSION} REQUIRED)

set(mrd_generated_files
    yardl/yardl.h
    protocols.cc
    types.cc
    binary/protocols.cc
)

set(mrd_model_dir ${mrd2_SOURCE_DIR}/model)
set(mrd_include_dir ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(mrd_generated_dir ${mrd_include_dir}/mrd)

file(GLOB mrd_model_files ${mrd_model_dir}/*.yml)
message(STATUS "mrd_model_files: ${mrd_model_files}")

list(TRANSFORM mrd_generated_files PREPEND ${mrd_generated_dir}/)
message(STATUS "mrd_generated_files: ${mrd_generated_files}")

add_custom_command(
    OUTPUT ${mrd_generated_files}
    COMMAND yardl generate -c cpp.overrideArrayHeader=gadgetron-arrays.h -c cpp.sourcesOutputDir=${mrd_generated_dir}
    WORKING_DIRECTORY ${mrd_model_dir}
    DEPENDS ${mrd_model_files}
    VERBATIM
)

add_library(pingvin_mrd SHARED ${mrd_generated_files} include/gadgetron-arrays.h)
add_library(mrd::mrd ALIAS pingvin_mrd)
set_target_properties(pingvin_mrd PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})
set_target_properties(pingvin_mrd PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(pingvin_mrd
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>"
    PUBLIC "$<BUILD_INTERFACE:${mrd_include_dir}>"
)

target_link_libraries(pingvin_mrd
    PUBLIC xtensor
    PUBLIC date::date
    PUBLIC nlohmann_json::nlohmann_json
    pingvin_toolbox_cpucore
)

install(TARGETS pingvin_mrd
    EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(DIRECTORY ${mrd_generated_dir}
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH}
    COMPONENT Devel
    FILES_MATCHING PATTERN "*.h"
)