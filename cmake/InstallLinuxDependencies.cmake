
#install dependencies
if (MKL_FOUND)
    if (HAS_64_BIT)
        set(MKL_REDIST_DIR ${MKLROOT_PATH}/mkl/lib/intel64)
        set(MKL_COMPILER_REDIST_DIR ${MKLROOT_PATH}/lib/intel64)
    else (HAS_64_BIT)
        set(MKL_REDIST_DIR ${MKLROOT_PATH}/mkl/lib/ia32)
        set(MKL_COMPILER_REDIST_DIR ${MKLROOT_PATH}/lib/ia32)
    endif (HAS_64_BIT)

    message("Install mkl libraries from ${MKL_REDIST_DIR} ")
    FILE(GLOB MKL_DLL ${MKL_REDIST_DIR}/*.so)
    foreach(fileName ${MKL_DLL})
        message("Install ${fileName} ")
        install( FILES ${fileName} DESTINATION lib COMPONENT main)
    endforeach(fileName)

    FILE(GLOB MKL_COMPILER_DLL ${MKL_COMPILER_REDIST_DIR}/libiomp5*.so)
    foreach(fileName ${MKL_COMPILER_DLL})
        message("Install ${fileName} ")
        install( FILES ${fileName} DESTINATION lib COMPONENT main)
    endforeach(fileName)
endif (MKL_FOUND)
