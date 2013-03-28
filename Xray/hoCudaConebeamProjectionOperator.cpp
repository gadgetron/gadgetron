#include "hoCudaConebeamProjectionOperator.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

#include "conebeam_projection.h"
#include <vector>

#include <stdio.h>
#include "hoCuNDArray_utils.h"
#include "hoNDArray_fileio.h"

//#define TIME_ITR
using namespace Gadgetron;
template<class REAL> void hoCudaConebeamProjectionOperator<REAL>
::mult_M( hoCuNDArray<REAL>* in1, hoCuNDArray<REAL>* out, bool accumulate )
{
  if( !preprocessed ){
  	BOOST_THROW_EXCEPTION(runtime_error("Error: hoCudaConebeamProjectionOperator : setup not performed"));
  }

  hoCuNDArray<REAL>* in = in1;
  hoCuNDArray<REAL> in2;
  /*  if (is_weights != NULL) {
    in2 = *in1;
    in  = &in2;
    *is_weights *= * in;
    }*/

  //unsigned int offset = 0;
  for (unsigned int b=0; b < binning->getBinningData().size(); b++) {
    //unsigned int numProjsInBin = projections_per_bin[b];
    /*
    //std::cout << "forward: " << b << std::endl;
    std::vector<unsigned int> ps_dims;
    ps_dims.push_back(ps_dims_in_pixels.x);
    ps_dims.push_back(ps_dims_in_pixels.y);
    ps_dims.push_back(numProjsInBin);
    hoCuNDArray<REAL> temp;
    unsigned int sliceSize = ps_dims_in_pixels.x * ps_dims_in_pixels.y;
    temp.create( &ps_dims, out->get_data_ptr() + offset*sliceSize );
    */

    floatd2 ps_dims_in_pixels_float = floatd2( ps_dims_in_pixels[0], ps_dims_in_pixels[1] );
    floatd2 ps_spacing_in_mm = geometry->getSpacingArray();
    floatd2 ps_dims_in_mm = ps_spacing_in_mm * ps_dims_in_pixels_float;

    floatd3 SAGx = geometry->getSAGxArray();
    floatd3 SAGy = geometry->getSAGyArray();
    float SDD = geometry->getSDD();
    float SAD = geometry->getSAD();
    conebeam_forwards_projection(*out, *in, b,
				      binning->getBinningData()[b],
				      angles, 
				      ppb, is_spacing_in_mm, ps_dims_in_mm, 
				      SAGx, SAGy, SDD, SAD, 
				      numSamplesPerRay, use_circular_cutoff, accumulate);

    //offset += projections_per_bin[b];
  }

  if (ps_weights != NULL) {
    *out *= * ps_weights;
  }

  if (ps_scale != 1.0) {
  	*out *= ps_scale;
  }
}

template<class REAL> void hoCudaConebeamProjectionOperator<REAL>
::mult_MH( hoCuNDArray<REAL>* in1, hoCuNDArray<REAL>* out, bool accumulate )
{
  if( !preprocessed ){
    BOOST_THROW_EXCEPTION(runtime_error("Error: hoCudaConebeamProjectionOperator : setup not performed"));
  }

  hoCuNDArray<REAL>* in = in1;
  hoCuNDArray<REAL> in2;
  if (ps_weights != NULL) {
    in2 = *in1;
    in  = &in2;
    *in *= *ps_weights;
  }

  for (unsigned int b=0; b < binning->getBinningData().size(); b++) {
    /*
      unsigned int numProjsInBin = projections_per_bin[b];
      //std::cout << "backward: " << b << std::endl;
      std::vector<unsigned int> ps_dims;
      ps_dims.push_back(ps_dims_in_pixels.x);
      ps_dims.push_back(ps_dims_in_pixels.y);
      ps_dims.push_back(numProjsInBin);
      hoCuNDArray<REAL> temp;
      unsigned int sliceSize = ps_dims_in_pixels.x * ps_dims_in_pixels.y;
      temp.create( &ps_dims, in->get_data_ptr() + offset*sliceSize );
    */
        
    floatd2 ps_dims_in_pixels_float = floatd2( ps_dims_in_pixels[0], ps_dims_in_pixels[1] );
    floatd2 ps_spacing_in_mm = geometry->getSpacingArray();
    floatd2 ps_dims_in_mm = ps_spacing_in_mm * ps_dims_in_pixels_float;

    floatd3 SAGx = geometry->getSAGxArray();
    floatd3 SAGy = geometry->getSAGyArray();
    float SDD = geometry->getSDD();
    float SAD = geometry->getSAD();
    conebeam_backwards_projection(*in, *out, b,
				       binning->getBinningData()[b],
				       angles, ppb, is_spacing_in_mm, ps_dims_in_mm,
				       SAGx, SAGy, SDD, SAD, 
				       use_circular_cutoff, use_fbp, accumulate);

  }

  if (is_weights != NULL) {
    *out *= *is_weights ;
  }

  if (is_scale != 1.0) {
  	*out *= is_scale;
  }

}

template<class REAL> void hoCudaConebeamProjectionOperator<REAL>
::mult_MH_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate )
{
  if( !preprocessed ){
  	BOOST_THROW_EXCEPTION(runtime_error("Error: hoCudaConebeamProjectionOperator : setup not performed"));
  }

  //if (gt!=NULL) delete gt;
  //gt = new GPUTimer("cgSolver Iteration");

  //GPUTimer *timer;
  //timer = new GPUTimer("hoCudaConebeamProjectionOperator mult_MH_M");

  // calc number of projections
  unsigned int numOfAllProjections = binning->getNumberOfProjections();

  /*
    unsigned int offset = 0;
    for (unsigned int b=0; b < projections_per_bin.size(); b++) {
    offset += projections_per_bin[b];
    }
  */

  // Make copy of input to avoid it being overwritten
  std::vector<unsigned int> ps_dims;
  ps_dims.push_back(ps_dims_in_pixels[0]);
  ps_dims.push_back(ps_dims_in_pixels[1]);
  ps_dims.push_back(numOfAllProjections);

  hoCuNDArray<REAL> temp;
  temp.create( &ps_dims );
        
#ifdef TIME_ITR
  GPUTimer *timer2;
  timer2 = new GPUTimer("forwards projection");
#endif

  mult_M(in, &temp, accumulate);
#ifdef TIME_ITR
  delete timer2;
  timer2 = new GPUTimer("backwards projection");
#endif

  mult_MH(&temp, out, accumulate);

#ifdef TIME_ITR
  delete timer2;
#endif
    

}

// Instantiations
template class hoCudaConebeamProjectionOperator<float>;
