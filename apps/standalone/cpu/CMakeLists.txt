include_directories( 
  ${Boost_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
)

#add_subdirectory(MRI)
add_subdirectory(denoising)
#add_subdirectory(deblurring)
