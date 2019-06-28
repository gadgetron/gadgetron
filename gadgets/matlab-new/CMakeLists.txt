
set(MATLAB_EXAMPLE_FILES ../matlab-new/examples/simple_recon.m)

message("Matlab example files: ${MATLAB_EXAMPLE_FILES}")

install(FILES ${MATLAB_EXAMPLE_FILES} DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
