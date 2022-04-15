# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/apps/standalone/cpu/registration/2d

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/apps/standalone/cpu/registration/2d/register_HS_2d_cpu")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/plplot:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/moco:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cartesian:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/dicom:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/grappa:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/spiral:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_noncartesian:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/gpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/hyper:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/fatwater:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/epi:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cmr:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/interventional_mri:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/plplot:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/python:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/pattern_recognition:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/ffd:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/nfft:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/nfft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/deblur:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fatwater:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_image:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/cmr:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/epi:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/spiral:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/solvers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/solvers/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/distributed:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/parallel:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/readers:/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_HS_2d_cpu")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/apps/standalone/cpu/registration/2d/register_CK_2d_cpu")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/plplot:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/moco:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cartesian:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/dicom:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/grappa:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/spiral:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_noncartesian:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/gpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/hyper:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/fatwater:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/epi:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cmr:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/interventional_mri:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/plplot:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/python:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/pattern_recognition:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/ffd:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/nfft:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/nfft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/deblur:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fatwater:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_image:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/cmr:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/epi:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/spiral:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/solvers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/solvers/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/distributed:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/parallel:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/readers:/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/register_CK_2d_cpu")
    endif()
  endif()
endif()

