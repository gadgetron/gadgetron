#include "parameterparser.h"
#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCudaConebeamProjectionOperator.h"
#include "hoNDArray_fileio.h"
#include "vector_td_utilities.h"
#include "GPUTimer.h"

#include <iostream>
#include <algorithm>
#include <sstream>

using namespace Gadgetron;
using namespace std;

int main(int argc, char** argv) 
{ 
  // Parse command line
  //

  ParameterParser parms(1024);
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input volume filename (.real)", true );
  parms.add_parameter( 'b', COMMAND_LINE_STRING, 1, "Binning filename (.hdf5) - 4D only", false );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output projections filename (.real)", true, "projections_simulated.real" );
  parms.add_parameter( 'h', COMMAND_LINE_STRING, 1, "Output acquisition filename (.h5)", true, "acquisition_simulated.h5" );
  parms.add_parameter( 'f', COMMAND_LINE_FLOAT, 3, "Input volume FOV in mm (3d)", true, "448, 448, 252" );
  parms.add_parameter( 'p', COMMAND_LINE_FLOAT, 2, "Projection plate size in pixels (2d)", true, "512, 256" );
  parms.add_parameter( 'q', COMMAND_LINE_FLOAT, 2, "Projection plate FOV in mm (2d)", true, "800.0, 400.0" );
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT, 1, "SAD", true, "1000.0" );
  parms.add_parameter( 's', COMMAND_LINE_FLOAT, 1, "SDD", true, "1500.0" );
  parms.add_parameter( 'u', COMMAND_LINE_FLOAT, 1, "Initial angle", true, "0.0" );
  parms.add_parameter( 'v', COMMAND_LINE_FLOAT, 1, "Angular spacing", true, "0.5" );
  parms.add_parameter( 'w', COMMAND_LINE_INT, 1, "Number of projections", true, "720" );
  parms.add_parameter( 'P', COMMAND_LINE_INT, 1, "Projections per batch", true, "50" );
  parms.add_parameter( 'S', COMMAND_LINE_INT, 1, "Samples per line integral", true, "0" );
  
  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ) {
    parms.print_parameter_list();
  }
  else{
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  std::string image_filename = (char*)parms.get_parameter('d')->get_string_value();
  std::string binning_filname = (char*)parms.get_parameter('b')->get_string_value();
  std::string projections_filename = (char*)parms.get_parameter('r')->get_string_value();
  std::string acquisition_filename = (char*)parms.get_parameter('h')->get_string_value();
  
  // Load volume
  //
  
  boost::shared_ptr< hoNDArray<float> > image = read_nd_array<float>( image_filename.c_str() );
  
  if( image->get_number_of_dimensions() < 3 ){
    std::cout << "Input image volume should have at least three dimensions" << std::endl;
    exit(1);
  }

  // Configuring...
  //

  uintd2 ps_dims_in_pixels( parms.get_parameter('p')->get_float_value(0), 
			    parms.get_parameter('p')->get_float_value(1) );
  
  floatd2 ps_dims_in_mm( parms.get_parameter('q')->get_float_value(0),
			 parms.get_parameter('q')->get_float_value(1) );

  floatd2 ps_spacing_in_mm = ps_dims_in_mm / ps_dims_in_pixels;

  float SAD = parms.get_parameter('a')->get_float_value();
  float SDD = parms.get_parameter('s')->get_float_value();

  uintd3 is_dims_in_pixels ( image->get_size(0),
			     image->get_size(1),
			     image->get_size(2) );
  
  floatd3 is_dims_in_mm( parms.get_parameter('f')->get_float_value(0), 
			 parms.get_parameter('f')->get_float_value(1), 
			 parms.get_parameter('f')->get_float_value(2) );

  floatd3 is_spacing_in_mm = is_dims_in_mm / is_dims_in_pixels;

  float start_angle = parms.get_parameter('u')->get_float_value();
  float angular_spacing = parms.get_parameter('v')->get_float_value();

  unsigned int number_of_projections = parms.get_parameter('w')->get_int_value();
  unsigned int projections_per_batch = parms.get_parameter('P')->get_int_value();

  unsigned int samples_per_ray = parms.get_parameter('S')->get_int_value();
  if( samples_per_ray == 0 ){
    float samples_per_pixel = 1.5f;
    samples_per_ray = (unsigned int)(samples_per_pixel*float(is_dims_in_mm[0]));
  }

  // Load the binning data (if provided)
  //
  
  boost::shared_ptr<CBCT_binning> binning( new CBCT_binning() );
  binning->set_as_default_3d_bin(number_of_projections);
  /*  if( binning_filename.size() > 0 ){
    binning = boost::shared_ptr<CBCT_binning>(new CBCT_binning());
    binning->load(bining_filename);
    } */
  binning->print();
  

  // Create projection angles array
  // - and offsets (assuming zeros for now)
  
  std::vector<float> angles;
  std::vector<floatd2> offsets;

  for( unsigned int i=0; i<number_of_projections; i++ ){
    float angle = start_angle + i*angular_spacing;
    angles.push_back(angle);
    offsets.push_back(floatd2(0.0f));
  }
  
  // Allocate array to hold the result
  //

  std::vector<unsigned int> ps_dims;
  ps_dims.push_back(ps_dims_in_pixels[0]);
  ps_dims.push_back(ps_dims_in_pixels[1]);
  ps_dims.push_back(number_of_projections);

  boost::shared_ptr< hoNDArray<float> > projections( new hoNDArray<float>(&ps_dims) );

  // Create geometry setup
  //

  boost::shared_ptr<CBCT_geometry> geometry( new CBCT_geometry() );
  geometry->set_SAD(SAD);
  geometry->set_SDD(SDD);
  geometry->set_spacing(ps_spacing_in_mm);
  geometry->set_angles(angles);
  geometry->set_offsets(offsets);

  // Create acquisition setup
  //

  boost::shared_ptr<CBCT_acquisition> acquisition( new CBCT_acquisition() );
  acquisition->set_geometry(geometry);
  acquisition->set_projections(projections);

  //
  // Standard 3D FDK reconstruction
  //
  
  boost::shared_ptr< hoCudaConebeamProjectionOperator > E( new hoCudaConebeamProjectionOperator() );
  E->setup( acquisition, binning, projections_per_batch, samples_per_ray, is_spacing_in_mm );

  {
    GPUTimer timer("Running CBCT forwards projection");
    E->mult_M( image.get(), projections.get() );
  }

  write_nd_array<float>( projections.get(), projections_filename.c_str() );
  acquisition->save( acquisition_filename );

  return 0;
}
