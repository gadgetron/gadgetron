add_subdirectory(pmri)

if (ARMADILLO_FOUND)
   add_subdirectory(epi)
   add_subdirectory(hyper)
endif (ARMADILLO_FOUND)

