if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()

add_subdirectory(cpu)

install(FILES
    KaiserBessel_kernel.h
    NFFT.h
    NFFT.hpp
    NFFTOperator.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
