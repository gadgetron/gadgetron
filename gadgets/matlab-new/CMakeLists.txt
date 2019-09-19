
set(MATLAB_EXAMPLE_FILES
        examples/simple_recon.m
        examples/bucket_recon.m)

message("Matlab example files: ${MATLAB_EXAMPLE_FILES}")

install(FILES ${MATLAB_EXAMPLE_FILES} DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
