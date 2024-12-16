add_subdirectory(cpu)

add_library(pingvin_toolbox_nfft INTERFACE)
target_include_directories(pingvin_toolbox_nfft
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES
    ConvolutionKernel.h
    ConvolutionKernel.hpp
    GriddingConvolution.h
    GriddingConvolution.hpp
    NFFT.h
    NFFT.hpp
    NFFTOperator.h
    NFFTOperator.hpp
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)