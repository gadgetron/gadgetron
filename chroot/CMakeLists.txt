configure_file("start-gadgetron.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron.sh @ONLY)
configure_file("enter-chroot-env.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/enter-chroot-env.sh @ONLY)
configure_file("gadgetron-dependency-query.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron-dependency-query.sh @ONLY)
configure_file("siemens_to_ismrmrd.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/siemens_to_ismrmrd.sh @ONLY)
configure_file("gadgetron_ismrmrd_client.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_ismrmrd_client.sh @ONLY)
configure_file("gt_alive.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/gt_alive.sh @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web_app.cfg @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web.conf.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web.conf @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_ld.conf.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web_ld.conf @ONLY)
configure_file("run-webapp.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-webapp.sh @ONLY)
configure_file("copy-cuda-lib.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/copy-cuda-lib.sh @ONLY)
configure_file("start-env.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-env.sh @ONLY)
configure_file("start.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start.sh @ONLY)
configure_file("start-webapp.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-webapp.sh @ONLY)
configure_file("run-gadgetron-dependency-query.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-gadgetron-dependency-query.sh @ONLY)
configure_file("run-gadgetron_ismrmrd_client.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-gadgetron_ismrmrd_client.sh @ONLY)
configure_file("run-gt_alive.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-gt_alive.sh @ONLY)
configure_file("run-siemens_to_ismrmrd.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-siemens_to_ismrmrd.sh @ONLY)
configure_file("mount_image.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/mount_image.sh @ONLY)
configure_file("start-gadgetron-from-image.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron-from-image.sh @ONLY)

set (LIBRARY_PATHS ${ISMRMRD_LIB_DIR})
set (LIBRARY_PATHS ${LIBRARY_PATHS}:${MKL_LIB_DIR})

FIND_LIBRARY( CUDA_LIBRARY NAMES libcuda.so )

string(FIND "${CUDA_LIBRARIES}" "libcudart.so" LIBCUDART_LOC)
string(FIND "${CUDA_LIBRARIES}" "libcuda.so" LIBCUDA_LOC)

if (LIBCUDA_LOC EQUAL -1)
    get_filename_component(CUDART_LIB_DIR ${CUDA_LIBRARIES} PATH)
else (LIBCUDA_LOC EQUAL -1)
    if(LIBCUDART_LOC LESS LIBCUDA_LOC)
        string(SUBSTRING "${CUDA_LIBRARIES}" 0 ${LIBCUDART_LOC} CUDART_LIB_DIR)
    else(LIBCUDART_LOC LESS LIBCUDA_LOC)
        string(SUBSTRING "${CUDA_LIBRARIES}" ${LIBCUDA_LOC} -1 CUDART_LIB_TEMP)
        string(SUBSTRING "${CUDART_LIB_TEMP}" 10 -1 CUDART_LIB_TEMP2)
        get_filename_component(CUDART_LIB_DIR ${CUDART_LIB_TEMP2} PATH)
    endif(LIBCUDART_LOC LESS LIBCUDA_LOC)
endif (LIBCUDA_LOC EQUAL -1)

get_filename_component(CUDA_LIB_DIR ${CUDA_LIBRARY} PATH)
get_filename_component(CUDAFFT_LIB_DIR ${CUDA_CUFFT_LIBRARIES} PATH)
get_filename_component(CUDABLAS_LIB_DIR ${CUDA_CUBLAS_LIBRARIES} PATH)

set (LIBRARY_PATHS ${LIBRARY_PATHS}:${CUDART_LIB_DIR}:${CUDA_LIB_DIR}:${CUDAFFT_LIB_DIR}:${CUDABLAS_LIB_DIR})

configure_file("copy-cuda-lib.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/copy-cuda-lib.sh @ONLY)

add_custom_target(chroot)
add_custom_target(chroot-local)

# Try to find GIT_SHA1 by calling gadgetron_info                                                                                           
if(EXISTS "${CMAKE_INSTALL_PREFIX}/gadgetron/bin/gadgetron_info")
  execute_process(COMMAND gadgetron_info
                  WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/gadgetron/bin
                  OUTPUT_VARIABLE GADGETRON_INFO
                  ERROR_VARIABLE GADGETRON_INFO_STDERR)

  string(REGEX MATCH "[0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z][0-9a-z]" GIT_SHA1 ${GADGETRON_INFO})

  if(GIT_SHA1)
    set(GADGETRON_GIT_SHA1 ${GIT_SHA1})
  endif(GIT_SHA1)

# Try to find GIT_SHA1 by calling git rev-parse HEAD
else(EXISTS "${CMAKE_INSTALL_PREFIX}/gadgetron/bin/gadgetron_info")
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GADGETRON_GIT_SHA1
                    ERROR_VARIABLE GADGETRON_GIT_STDERR)

    string(STRIP "${GADGETRON_GIT_SHA1}" GADGETRON_GIT_SHA1)
    string(LENGTH "${GADGETRON_GIT_SHA1}" GADGETRON_GIT_SHA1_LEN)
    if(${GADGETRON_GIT_SHA1_LEN} LESS 40)
      message(WARNING "Could not determine SHA-1 hash: ${GADGETRON_GIT_STDERR}")
      set(GADGETRON_GIT_SHA1 "NA")
    endif(${GADGETRON_GIT_SHA1_LEN} LESS 40)
  else(GIT_FOUND)
    set(GADGETRON_GIT_SHA1 "NA")
  endif(GIT_FOUND)
endif(EXISTS "${CMAKE_INSTALL_PREFIX}/gadgetron/bin/gadgetron_info")

add_custom_command(
  TARGET chroot
  POST_BUILD
  VERBATIM
  COMMENT "Creating chroot tar and img file"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND sudo bash create_chroot.sh 
  ${CMAKE_INSTALL_PREFIX}
  ${CMAKE_BINARY_DIR}
  ${GADGETRON_GIT_SHA1}
  ${LIBRARY_PATHS}
  )

set(PACKAGES_PATH "/home/ubuntu/packages")

add_custom_command(
  TARGET chroot-local
  POST_BUILD
  VERBATIM
  COMMENT "Creating chroot tar and img file using local packages"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND sudo bash create_chroot-local.sh
  ${CMAKE_INSTALL_PREFIX}
  ${CMAKE_BINARY_DIR}
  ${GADGETRON_GIT_SHA1}
  ${LIBRARY_PATHS}
  ${PACKAGES_PATH}
  )

install(PROGRAMS 
  ${CMAKE_CURRENT_BINARY_DIR}/copy-cuda-lib.sh 
  ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron.sh
  ${CMAKE_CURRENT_BINARY_DIR}/enter-chroot-env.sh
  ${CMAKE_CURRENT_BINARY_DIR}/gadgetron-dependency-query.sh
  ${CMAKE_CURRENT_BINARY_DIR}/siemens_to_ismrmrd.sh
  ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_ismrmrd_client.sh
  ${CMAKE_CURRENT_BINARY_DIR}/gt_alive.sh
  ${CMAKE_CURRENT_BINARY_DIR}/run-webapp.sh
  ${CMAKE_CURRENT_BINARY_DIR}/run-gadgetron-dependency-query.sh
  ${CMAKE_CURRENT_BINARY_DIR}/run-gadgetron_ismrmrd_client.sh
  ${CMAKE_CURRENT_BINARY_DIR}/run-gt_alive.sh
  ${CMAKE_CURRENT_BINARY_DIR}/run-siemens_to_ismrmrd.sh
  ${CMAKE_CURRENT_BINARY_DIR}/start-env.sh
  ${CMAKE_CURRENT_BINARY_DIR}/start.sh
  ${CMAKE_CURRENT_BINARY_DIR}/start-webapp.sh
  ${CMAKE_CURRENT_BINARY_DIR}/mount_image.sh
  ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron-from-image.sh
  ${CMAKE_SOURCE_DIR}/chroot/mount.sh
  ${CMAKE_SOURCE_DIR}/chroot/stop.sh
  ${CMAKE_SOURCE_DIR}/chroot/umount_image.sh
  ${CMAKE_SOURCE_DIR}/chroot/install_chroot_image.sh
  DESTINATION ${GADGETRON_INSTALL_CHROOT_SCRIPTS_PATH}
  COMPONENT scripts
  )
