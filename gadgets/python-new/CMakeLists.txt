
set(PYTHON_EXAMPLE_FILES
        examples/simple_recon.py
        examples/simple_bucket_recon.py)

message("Python example files: ${PYTHON_EXAMPLE_FILES}")

install(FILES ${PYTHON_EXAMPLE_FILES} DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

