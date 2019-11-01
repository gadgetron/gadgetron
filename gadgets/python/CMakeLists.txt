
set(gadgetron_python_config_files
		legacy/config/pseudoreplica.xml
		legacy/config/python_buckets.xml
		legacy/config/python_short.xml
		legacy/config/python_image_array.xml
		legacy/config/python_image_array_recon.xml
		legacy/config/python_passthrough.xml
		legacy/config/python_short.xml
		config/Generic_Cartesian_Grappa_RealTimeCine_Python.xml)

set(gadgetron_python_gadgets_files
		legacy/gadgets/accumulate_and_recon.py
		legacy/gadgets/array_image.py
		legacy/gadgets/bucket_recon.py
        legacy/gadgets/image_array_recon.py
		legacy/gadgets/image_array_recon_rtcine_plotting.py
		legacy/gadgets/passthrough.py
		legacy/gadgets/passthrough_array_image.py
		legacy/gadgets/pseudoreplicagather.py
		legacy/gadgets/remove_2x_oversampling.py
		legacy/gadgets/rms_coil_combine.py)

source_group(config FILES ${gadgetron_python_config_files})
source_group(gadgets FILES ${gadgetron_python_gadgets_files})

install(FILES ${gadgetron_python_gadgets_files}
        DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(FILES ${gadgetron_python_config_files}
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
