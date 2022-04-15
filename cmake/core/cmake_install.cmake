# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/core

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/javeda2/anaconda3/envs/gadgetron")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/core/libgadgetron_core.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/core/libgadgetron_core.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/core/libgadgetron_core.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_core.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/core/Channel.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Channel.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/core/ChannelIterator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Message.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Message.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/core/MPMCChannel.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Gadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Context.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Gadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Reader.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Types.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Types.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/core/TypeTraits.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Writer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Writer.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Node.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/PureGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/LegacyACE.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/PropertyMixin.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/GadgetContainerMessage.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/ChannelAlgorithms.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/variant.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/core/MessageID.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Storage.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/Process.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron/io" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/core/io/adapt_struct.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/io/from_string.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/io/ismrmrd_types.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/io/primitives.h"
    "/home/javeda2/mrprogs/gadgetron_fork/core/io/primitives.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/config" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/core/config/distributed_default.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/core/config/distributed_generic_default.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/core/config/distributed_image_default.xml"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/core/readers/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/core/parallel/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/core/distributed/cmake_install.cmake")

endif()

