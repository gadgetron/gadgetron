#include "parameterparser.h"
#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCuConebeamProjectionOperator.h"
#include "hoNDArray_fileio.h"
#include "vector_td_utilities.h"
#include "GPUTimer.h"
#include "setup_grid.h"

#include <iostream>
#include <algorithm>
#include <sstream>

using namespace Gadgetron;
using namespace std;

// Utility to load offsets - if file is provided - from an HDF5 file
//

std::vector<floatd2> 
get_offsets( std::string filename )
{  
  hsize_t dim;
  
  hid_t file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
  
  herr_t errCode = H5LTget_dataset_info(file_id,"/offsetx",&dim,NULL,NULL);
  if (errCode < 0) 
    throw std::runtime_error("Error getting /offsetx dataset info from file.");

  std::vector<float> offsets_x = std::vector<float>(dim,0.0f);
  errCode=H5LTread_dataset (file_id, "/offsetx", H5T_NATIVE_FLOAT, &offsets_x[0]);
  if (errCode < 0)
    throw std::runtime_error("Error reading /offsetx from file.");
  
  errCode=H5LTget_dataset_info(file_id,"/offsety",&dim,NULL,NULL);
  if (errCode < 0)
    throw std::runtime_error("Error getting /offsety dataset info from file.");

  std::vector<float> offsets_y = std::vector<float>(dim,0.0f);

  errCode = H5LTread_dataset (file_id, "/offsety", H5T_NATIVE_FLOAT, &offsets_y[0]);
  if (errCode < 0)
    throw std::runtime_error("Error reading /offsety from file.");
  
  if( offsets_x.size() != offsets_y.size() ){
    throw std::runtime_error("CBCT_geometry::load : x/y offset arrays has different lengths");
  }

  std::vector<floatd2> res;
  for( unsigned int i=0; i<offsets_x.size(); i++ )
    res.push_back(floatd2( offsets_x[i], offsets_y[i]));
  
  return res;
}

int main(int argc, char** argv) 
{ 
  // Parse command line
  //

  ParameterParser parms(1024);
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input volume filename (.real)", true );
  parms.add_parameter( 'b', COMMAND_LINE_STRING, 1, "Binning filename (.h5)", false );
  parms.add_parameter( 'o', COMMAND_LINE_STRING, 1, "Offsets filename (.h5)", false );
  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output projections filename (.real)", true, "projections_simulated.real" );
  parms.add_parameter( 'h', COMMAND_LINE_STRING, 1, "Output acquisition filename (.h5)", true, "acquisition_simulated.h5" );
  parms.add_parameter( 'f', COMMAND_LINE_FLOAT, 3, "Input volume FOV in mm (3d)", true, "448, 448, 252" );
  parms.add_parameter( 'p', COMMAND_LINE_FLOAT, 2, "Projection plate size in pixels (2d)", true, "512, 256" );
  parms.add_parameter( 'q', COMMAND_LINE_FLOAT, 2, "Projection plate FOV in mm (2d)", true, "800.0, 400.0" );
  parms.add_parameter( 'a', COMMAND_LINE_FLOAT, 1, "SAD", true, "1000.0" );
  parms.add_parameter( 's', COMMAND_LINE_FLOAT, 1, "SDD", true, "1500.0" );
  parms.add_parameter( 'u', COMMAND_LINE_FLOAT, 1, "Initial angle (degrees)", true, "0.0" );
  parms.add_parameter( 'v', COMMAND_LINE_FLOAT, 1, "Angular spacing (degrees)", true, "0.5" );
  parms.add_parameter( 'w', COMMAND_LINE_INT, 1, "Number of projections", true, "720" );
  parms.add_parameter( 'P', COMMAND_LINE_INT, 1, "Projections per batch", false );
  parms.add_parameter( 'S', COMMAND_LINE_FLOAT, 1, "Samples per pixel (float) in integral", false );
  
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
  std::string projections_filename = (char*)parms.get_parameter('r')->get_string_value();
  std::string acquisition_filename = (char*)parms.get_parameter('h')->get_string_value();
  
  // Load volume
  //
  
  boost::shared_ptr< hoCuNDArray<float> > image(new hoCuNDArray<float>(*read_nd_array<float>( image_filename.c_str() )));
  
  if( image->get_number_of_dimensions() < 3 ){
    std::cout << "Input image volume should have at least three dimensions" << std::endl;
    exit(1);
  }

  // Add default temporal dimension of 1 since the operator only takes four-dimensional images
  //

  if( image->get_number_of_dimensions() == 3 ){
    std::vector<unsigned int> dims = *image->get_dimensions();
    dims.push_back(1);
    image->reshape(&dims);
  }


  // Configuring...
  //

  uintd2 ps_dims_in_pixels( parms.get_parameter('p')->get_float_value(0), 
			    parms.get_parameter('p')->get_float_value(1) );
  
  floatd2 ps_dims_in_mm( parms.get_parameter('q')->get_float_value(0),
			 parms.get_parameter('q')->get_float_value(1) );

  float SAD = parms.get_parameter('a')->get_float_value();
  float SDD = parms.get_parameter('s')->get_float_value();

  uintd3 is_dims_in_pixels ( image->get_size(0),
			     image->get_size(1),
			     image->get_size(2) );
  
  floatd3 is_dims_in_mm( parms.get_parameter('f')->get_float_value(0), 
			 parms.get_parameter('f')->get_float_value(1), 
			 parms.get_parameter('f')->get_float_value(2) );

  float start_angle = parms.get_parameter('u')->get_float_value();
  float angular_spacing = parms.get_parameter('v')->get_float_value();

  unsigned int number_of_projections = parms.get_parameter('w')->get_int_value();

  // Load or generate binning data
  //
  
  boost::shared_ptr<CBCT_binning> binning( new CBCT_binning() );

  if (parms.get_parameter('b')->get_is_set()){
    std::string binningdata_filename = (char*)parms.get_parameter('b')->get_string_value();
    std::cout << "Using binning data file: " << binningdata_filename << std::endl;
    binning->load(binningdata_filename);

    if( binning->get_maximum_projection_index() >= number_of_projections ) {
      std::cout << "Maximum projection index in binning file (" << 
	binning->get_maximum_projection_index() << 
	") exceeds the number of projections requested at the command line (" << 
	number_of_projections <<
	")" << std::endl;
      exit(1);
    }
  } 
  else 
    binning->set_as_default_3d_bin(number_of_projections);

  binning->print();
  
  // Create projection angles array
  //
  
  std::vector<float> angles;

  for( unsigned int i=0; i<number_of_projections; i++ ){
    float angle = start_angle + i*angular_spacing;
    angles.push_back(angle);
  }
  
  // Create projection offsets array
  //

  std::vector<floatd2> offsets;

  if (parms.get_parameter('o')->get_is_set()){
    std::string offsets_filename = (char*)parms.get_parameter('o')->get_string_value();
    std::cout << "Using offsets file: " << offsets_filename << std::endl;
    offsets = get_offsets(offsets_filename);
  } 
  else{
    for( unsigned int i=0; i<number_of_projections; i++ ){
      offsets.push_back(floatd2(0.0f));
    } 
  }   
  
  // Allocate and clear array to hold the result
  //
  
  std::vector<unsigned int> ps_dims;
  ps_dims.push_back(ps_dims_in_pixels[0]);
  ps_dims.push_back(ps_dims_in_pixels[1]);
  ps_dims.push_back(number_of_projections);

  boost::shared_ptr< hoCuNDArray<float> > projections( new hoCuNDArray<float>(&ps_dims) );
  clear(projections.get()); // Since the binning might not write to all projections

  // Create geometry setup
  //

  boost::shared_ptr<CBCT_geometry> geometry( new CBCT_geometry() );
  geometry->set_SAD(SAD);
  geometry->set_SDD(SDD);
  geometry->set_FOV(ps_dims_in_mm);
  geometry->set_angles(angles);
  geometry->set_offsets(offsets);

  // Create acquisition setup
  //

  boost::shared_ptr<CBCT_acquisition> acquisition( new CBCT_acquisition() );
  acquisition->set_geometry(geometry);
  acquisition->set_projections(projections);

  // Define conebeam projection operator
  // - and configure based on input parameters
  //

  boost::shared_ptr< hoCuConebeamProjectionOperator > E( new hoCuConebeamProjectionOperator() );
  
  CommandLineParameter *parm = parms.get_parameter('P');
  if( parm && parm->get_is_set() )
    E->set_num_projections_per_batch( parm->get_int_value() );
  
  parm = parms.get_parameter('S');  
  if( parm && parm->get_is_set() ) 
    E->set_num_samples_per_pixel( parm->get_float_value() );
  
  E->setup( acquisition, binning, is_dims_in_mm );

  // Initialize the device
  // - just to report more accurate timings
  //

  cudaThreadSynchronize();

  //
  // Forwards projection (X-ray image simulation)
  //
  
  {
    GPUTimer timer("Running CBCT forwards projection");
    E->mult_M( image.get(), projections.get() );
    cudaThreadSynchronize();
  }

  write_nd_array<float>( projections.get(), projections_filename.c_str() );
  acquisition->save( acquisition_filename );

  return 0;
}
