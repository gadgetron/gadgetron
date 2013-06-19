#pragma once

#include "hoCuNDArray.h"

/*
  Forwards and backwards projection
  ---------------------------------
  RealFloatArray *projections: four-dimensional ndarray (2d projection x num_coils x num_projections)
  RealFloatArray       *image: four-dimensional ndarray (3d volume x num_coils)
*/

namespace Gadgetron{

// Forwards projection of the 3D volume of a given phase onto projections
template <class TYPE>
bool conebeam_forwards_projection( hoCuNDArray<TYPE>& projections, hoCuNDArray<TYPE>& image, unsigned int bin,
                                   std::vector<unsigned int> binningdata,
                                   std::vector<float> angles,
                                   unsigned int ppb, floatd3 is_box_dims, floatd2 ip_dims,
                                   std::vector<float> offsetx, std::vector<float> offsety,
                                   float SDD, float SAD,
                                   unsigned int numSamplesPerRay,  bool accumulate);

// Backproject all projections for a given phase
template <class TYPE>
bool conebeam_backwards_projection( hoCuNDArray<TYPE>& projections, hoCuNDArray<TYPE>& image, unsigned int bin,
									std::vector<unsigned int> binningdata,
                                    std::vector<float> angles,
                                    unsigned int ppb, floatd3 is_box_dims, floatd2 ip_dims,
                                    std::vector<float> offsetx, std::vector<float> offsety,
                                    float SDD, float SAD,
                                    bool use_fbp, bool accumulate);
}
