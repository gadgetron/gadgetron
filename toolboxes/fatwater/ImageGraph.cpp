#include "ImageGraph.h"

template<> Gadgetron::vector_td<int, 2> Gadgetron::ImageGraph<2>::index_to_offset[4] = {Gadgetron::vector_td<int,2>(-1,0),
                                                                             Gadgetron::vector_td<int,2>(1,0),
                                                                             Gadgetron::vector_td<int,2>(0,-1),
                                                                             Gadgetron::vector_td<int,2>(0,1)};


template<> Gadgetron::vector_td<int, 3> Gadgetron::ImageGraph<3>::index_to_offset[6] = {Gadgetron::vector_td<int,3>(-1,0,0),
                                                                             Gadgetron::vector_td<int,3>(1,0,0),
                                                                             Gadgetron::vector_td<int,3>(0,-1,0),
                                                                             Gadgetron::vector_td<int,3>(0,1,0),
                                                                             Gadgetron::vector_td<int,3>(0,0,-1),
                                                                             Gadgetron::vector_td<int,3>(0,0,1)

};


