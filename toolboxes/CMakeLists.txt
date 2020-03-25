if (WIN32)
    add_definitions(-DTIXML_USE_STL)
endif ()


add_subdirectory(log)

add_subdirectory(operators)
add_subdirectory(solvers)

add_subdirectory(fft)
add_subdirectory(core)
add_subdirectory(mri_core)
add_subdirectory(mri)
add_subdirectory(cmr)

add_subdirectory(mri_image)

add_subdirectory(klt)
add_subdirectory(fatwater)
add_subdirectory(deblur)
add_subdirectory(nfft)
add_subdirectory(dwt)
add_subdirectory(registration)

add_subdirectory(ffd)
add_subdirectory(image)
add_subdirectory(pattern_recognition)
add_subdirectory(denoise)

add_subdirectory(image_io)

# Should we compile the python toolbox
if (BUILD_PYTHON_SUPPORT)
    message("Find python; build gadgetron toolbox python")
    add_subdirectory(python)
else ()
    message("Cannot find python; do not build gadgetron toolbox python")
endif ()

if (ARMADILLO_FOUND AND Boost_FOUND AND PLPLOT_FOUND)
    add_subdirectory(plplot)
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)