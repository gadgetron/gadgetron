find_package(Ismrmrd REQUIRED)
find_package(Octave REQUIRED)

include_directories(   ${OCTAVE_INCLUDE_DIRS}
					   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR} 
				       ${ISMRMRD_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

link_directories(${OCTAVE_LINK_DIRS})

add_library(gadgetronoctavecommunicator SHARED OctaveCommunicator.cpp)
target_link_libraries(gadgetronoctavecommunicator ${OCTAVE_LIBRARY} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

add_library(GadgetronReturnIsmrmrdAcquisition MODULE GadgetronReturnIsmrmrdAcquisition.cpp)
SET_TARGET_PROPERTIES(GadgetronReturnIsmrmrdAcquisition PROPERTIES SUFFIX .oct PREFIX "")
target_link_libraries(GadgetronReturnIsmrmrdAcquisition ${OCTAVE_LIBRARY} gadgetronoctavecommunicator optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

add_library(XMLGetXPath MODULE pugixml.cpp XMLGetXPath.cpp)
SET_TARGET_PROPERTIES(XMLGetXPath PROPERTIES SUFFIX .oct PREFIX "")
target_link_libraries(XMLGetXPath ${OCTAVE_LIBRARY})

add_library(GadgetronReturnIsmrmrdImage MODULE GadgetronReturnIsmrmrdImage.cpp)
SET_TARGET_PROPERTIES(GadgetronReturnIsmrmrdImage PROPERTIES SUFFIX .oct PREFIX "")
target_link_libraries(GadgetronReturnIsmrmrdImage ${OCTAVE_LIBRARY} gadgetronoctavecommunicator optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

add_library(gadgetronoctave SHARED OctaveGadget.cpp)
target_link_libraries(gadgetronoctave ${OCTAVE_LIBRARY} gadgetronoctavecommunicator cpucore)
	
install(TARGETS gadgetronoctavecommunicator DESTINATION lib)
install(TARGETS gadgetronoctave DESTINATION lib)
install(TARGETS GadgetronReturnIsmrmrdAcquisition DESTINATION octave)
install(TARGETS GadgetronReturnIsmrmrdImage DESTINATION octave)
install(TARGETS XMLGetXPath DESTINATION octave)

install(FILES octave/my_gadget_reference.m 
              octave/my_recon_function.m 
              octave/my_config_function.m 
              octave/gadget_reference_downsample_2x.m 
              octave/downsample_2x.m 
              octave/configure_downsample_2x.m 
              octave/ismrm_transform_kspace_to_image.m 
              octave/ismrm_transform_image_to_kspace.m 
              octave/gadget_reference_accumulator.m 
              octave/accumulator.m 
              octave/configure_accumulator.m 
        DESTINATION octave)

install(FILES octave.xml DESTINATION config)
