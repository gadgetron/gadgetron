if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_PYTHON__)
endif ()

# TODO: The Gadgetron Python code uses Numpy C-API code

set(gadgetron_python_header_files
        GadgetReference.h
        gadgetronpython_export.h
        PythonGadget.h)

set(gadgetron_python_src_files GadgetronPythonMRI.cpp
        GadgetReference.cpp
        PythonGadget.cpp)

set(gadgetron_python_config_files
        config/pseudoreplica.xml
        config/python.xml
        config/python_buckets.xml
        config/python_ideal_cg.xml
        config/python_short.xml
        config/python_image_array.xml
        config/python_image_array_recon.xml
        config/python_passthrough.xml
        config/python_short.xml
        config/python_tpat_snr_scale.xml
        config/Generic_Cartesian_Grappa_RealTimeCine_Python.xml
        )

set(gadgetron_python_utils_files
        utils/gadgetron_python_to_xml.py
        utils/gadgetron_run_python_chain.py
        utils/gadgetron_xml_to_python.py)

set(gadgetron_python_gadgets_files
        gadgets/accumulate_and_recon.py
        gadgets/bucket_recon.py
        gadgets/IDEAL.py
        gadgets/image_viewer.py
        gadgets/passthrough.py
        gadgets/pseudoreplicagather.py
        gadgets/remove_2x_oversampling.py
        gadgets/rms_coil_combine.py
        gadgets/tpat_snr_scale.py
        gadgets/array_image.py
        gadgets/image_array_recon.py
        gadgets/passthrough_array_image.py
        gadgets/image_array_recon_rtcine_plotting.py
        )

set(gadgetron_python_examples_files
        ../python-new/examples/simple_recon.py)

source_group(config FILES ${gadgetron_python_config_files})
source_group(utils FILES ${gadgetron_python_utils_files})
source_group(gadgets FILES ${gadgetron_python_gadgets_files})
source_group(examples FILES ${gadgetron_python_examples_files})
source_group(models             FILES      ${gadgetron_python_models_files})

add_library(GadgetronPythonMRI MODULE
        GadgetronPythonMRI.cpp
        GadgetReference.h
        GadgetReference.cpp
        )

#We should probably not set soversion on the python module. Causes problems with clang
#set_target_properties(GadgetronPythonMRI PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

add_library(gadgetron_python SHARED
        ${gadgetron_python_header_files}
        ${gadgetron_python_src_files}
        ${gadgetron_python_config_files}
        ${gadgetron_python_utils_files}
        ${gadgetron_python_gadgets_files}
        ${gadgetron_python_examples_files}${gadgetron_python_models_files} )

set_target_properties(gadgetron_python PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_python
        gadgetron_core
        gadgetron_toolbox_python
        gadgetron_toolbox_log)

target_include_directories(gadgetron_python
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(GadgetronPythonMRI
        gadgetron_python
        boost
        python )
if (WIN32)
    set_target_properties(GadgetronPythonMRI PROPERTIES SUFFIX .pyd)
    set_target_properties(gadgetron_python PROPERTIES LINK_FLAGS "/LIBPATH:${PYTHON_INCLUDE_DIR}/../libs")
endif ()

set_target_properties(GadgetronPythonMRI PROPERTIES PREFIX "")


install(TARGETS gadgetron_python
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)
install(TARGETS GadgetronPythonMRI DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(FILES
        ${gadgetron_python_gadgets_files}
        ${gadgetron_python_examples_files}
        ${gadgetron_python_models_files}
        DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(FILES ${gadgetron_python_config_files}
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

add_subdirectory(utils)
