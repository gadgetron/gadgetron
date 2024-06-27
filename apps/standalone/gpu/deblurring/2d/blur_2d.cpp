/*
  Example code to blur an image and generate input data for the deblurring apps.
*/

#include "hoNDArray_fileio.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "parameterparser.h"
#include "cuConvolutionOperator.h"

#include <iostream>
#include <math.h>

using namespace std;
using namespace Gadgetron;

// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;

int main( int argc, char** argv) 
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input image file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output blurred image file name (.cplx)", true, "blurred_image.cplx" );
  parms.add_parameter( 'k', COMMAND_LINE_STRING, 1, "Output kernel image file name (.cplx)", true, "kernel_image.cplx" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    GINFO_STREAM(" Running with the following parameters: " << endl);
    parms.print_parameter_list();
  }
  else{
    GINFO_STREAM(" Some required parameters are missing: " << endl);
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  // Load image from disk (single precision assumed)
  boost::shared_ptr< hoNDArray<float> > _host_image = 
    read_nd_array<float>((char*)parms.get_parameter('d')->get_string_value());

  if( !(_host_image->get_number_of_dimensions() == 2) ){
    GINFO_STREAM(endl << "Input image is not two-dimensional. Quitting.\n" << endl);
    return 1;
  }

  // Convert to _real
  hoNDArray<_real> host_image; host_image.create(_host_image->dimensions());
  for( unsigned int i=0; i<host_image.get_number_of_elements(); i++ )
    host_image.get_data_ptr()[i] = (_real) _host_image->get_data_ptr()[i];
    
  // Upload host image to device, normalize, and convert to complex type
  cuNDArray<_real> _image(host_image);
  normalize( &_image, _real(1) );
  boost::shared_ptr< cuNDArray<_complext> > image = real_to_complex<_complext>( &_image );
  
  // Setup resulting blurred image
  cuNDArray<_complext> blurred_image; 
  blurred_image.create(image->get_dimensions());
  
  // Generate convolution kernel (just do this on the host for now)
  _real sigma = 2.5;
  hoNDArray<_real> host_kernel;
  host_kernel.create(image->dimensions());
  for( unsigned int y=0; y<image->get_size(1); y++ ){
    for( unsigned int x=0; x<image->get_size(0); x++ ){
      _real biasx = (_real)(image->get_size(0)>>1);
      _real biasy = (_real)(image->get_size(1)>>1);
      _real cx = (_real)x-biasx;
      _real cy = (_real)y-biasy;
      host_kernel.get_data_ptr()[y*image->get_size(0)+x] = 1.0/(2.0*M_PI*sigma*sigma)*exp(-1.0*((cx*cx)/(2.0*sigma*sigma)+(cy*cy)/(2.0*sigma*sigma)));
    }
  }

  cuNDArray<_real> _kernel(host_kernel);
  boost::shared_ptr< cuNDArray<_complext> > kernel = real_to_complex<_complext>( &_kernel );

  // Normalize kernel
  _real scale = asum(kernel.get());
  *kernel /= scale;

  // Create convolution operator and assign kernel
  cuConvolutionOperator<_real,2> conv;
  conv.set_kernel( kernel.get() );  

  // Convolve
  conv.mult_M( image.get(), &blurred_image );

  //
  // Output result
  //
  
  boost::shared_ptr< hoNDArray<_complext> > blurred_image_host = blurred_image.to_host();
  write_nd_array<_complext>( blurred_image_host.get(), (char*)parms.get_parameter('r')->get_string_value());

  boost::shared_ptr< hoNDArray<_real> > host_norm = abs(&blurred_image)->to_host();
  write_nd_array<_real>( host_norm.get(), "blurred_image.real" );

  boost::shared_ptr< hoNDArray<_complext> > kernel_image_host = kernel->to_host();
  write_nd_array<_complext>( kernel_image_host.get(), (char*)parms.get_parameter('k')->get_string_value());

  return 0;
}
