if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()

add_subdirectory(cpu)

add_library(gadgetron_toolbox_nfft INTERFACE)
target_include_directories(gadgetron_toolbox_nfft
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES
    nfft_export.h
    KaiserBessel_kernel.h
    NFFT.h
    NFFT.hpp
    NFFTOperator.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
