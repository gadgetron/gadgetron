/*
  Example code to blur an image and generate input data for the deblurring apps.
*/

#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "parameterparser.h"
#include "cuConvolutionOperator.h"
#include "ndarray_vector_td_utilities.h"

#include <iostream>
#include <math.h>

using namespace std;

// Define desired precision (note that decent deblurring of noisy images requires double precision)
typedef float _real; 
typedef complext<_real>::Type _complext;

int main( int argc, char** argv) 
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input image file name (.real)", true );
  parms.add_parameter( 'k', COMMAND_LINE_STRING, 1, "In kernel image file name (.real)", true );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image file name (.cplx)", true, "blurred_image.cplx" );
  parms.add_parameter( 'K', COMMAND_LINE_STRING, 1, "Output kernel file name (.cplx)", true, "kernel_image.cplx" );

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  // Load image and kernel from disk (single precision assumed)
  //
  boost::shared_ptr< hoNDArray<float> > _host_image = 
    read_nd_array<float>((char*)parms.get_parameter('d')->get_string_value());

  boost::shared_ptr< hoNDArray<float> > _host_kernel = 
    read_nd_array<float>((char*)parms.get_parameter('k')->get_string_value());

  if( !(_host_image->get_number_of_dimensions() == 3) ){
    cout << endl << "Input image is not three-dimensional. Quitting.\n" << endl;
    return 1;
  }

  if( !(_host_kernel->get_number_of_dimensions() == 3) ){
    cout << endl << "Input kernel is not three-dimensional. Quitting.\n" << endl;
    return 1;
  }

  // Convert image and kernel to _real
  //
  hoNDArray<_real> host_image; host_image.create(_host_image->get_dimensions().get()); 
  for( unsigned int i=0; i<host_image.get_number_of_elements(); i++ )
    host_image.get_data_ptr()[i] = (_real) _host_image->get_data_ptr()[i];
    
  hoNDArray<_real> host_kernel; host_kernel.create(_host_kernel->get_dimensions().get()); 
  for( unsigned int i=0; i<host_kernel.get_number_of_elements(); i++ )
    host_kernel.get_data_ptr()[i] = (_real) _host_kernel->get_data_ptr()[i];

  // Upload host image and to device, normalize, and convert to complex type
  //
  cuNDArray<_real> _image(&host_image);
  cuNDA_normalize( &_image, get_one<_real>() );
  boost::shared_ptr< cuNDArray<_complext> > image = cuNDA_real_to_complext<_real>( &_image );
  
  cuNDArray<_real> _kernel(&host_kernel);
  cuNDA_normalize( &_kernel, get_one<_real>() );
  boost::shared_ptr< cuNDArray<_complext> > kernel = cuNDA_real_to_complext<_real>( &_kernel );

  // Setup resulting blurred image
  cuNDArray<_complext> blurred_image;
  blurred_image.create(image->get_dimensions().get());
  
  // Create convolution operator and assign kernel
  cuConvolutionOperator<_real> conv;
  conv.set_kernel( kernel.get() );  

  // Convolve
  conv.mult_M( image.get(), &blurred_image );

  //
  // Output result
  //
  
  boost::shared_ptr< hoNDArray<_complext> > blurred_image_host = blurred_image.to_host();
  write_nd_array<_complext>( blurred_image_host.get(), (char*)parms.get_parameter('r')->get_string_value());

  boost::shared_ptr< hoNDArray<_real> > host_norm = cuNDA_norm<_real>(&blurred_image)->to_host();
  write_nd_array<_real>( host_norm.get(), "blurred_image.real" );

  boost::shared_ptr< hoNDArray<_complext> > kernel_image_host = kernel->to_host();
  write_nd_array<_complext>( kernel_image_host.get(), (char*)parms.get_parameter('K')->get_string_value());

  boost::shared_ptr< hoNDArray<_real> > host_norm_kernel = cuNDA_norm<_real>(kernel.get())->to_host();
  write_nd_array<_real>( host_norm_kernel.get(), "kernel_image.real" );

  return 0;
}
