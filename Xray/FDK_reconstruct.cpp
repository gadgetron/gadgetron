#include "hoNDArray_utils.h"
#include "radial_utilities.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "parameterparser.h"

#include <iostream>
#include "math_constants.h"
#include <algorithm>

#include "imageOperator.h"
#include "identityOperator.h"

#include "hoCuPartialDerivativeOperator.h"

#include "hoCudaConebeamProjectionOperator.h"

// TMP (blurring)
#include "cuConvolutionOperator.h"

#include "hoCuNDArray_blas.h"
#include "hoCuNDArray_operators.h"
#include "hoNDArray_blas.h"
#include "cgSolver.h"
#include "PS_Dataset.h"
#include <sstream>
#include "complext.h"

#include "encodingOperatorContainer.h"

using namespace Gadgetron;
// Define desired precision
typedef float _real; 
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;





int main(int argc, char** argv) {


	std::string cmd = "";
	for (unsigned int i=0; i<argc; i++) {
			if (i!=0) cmd += " ";
			cmd += argv[i];
	}
  // Parse command line
  ParameterParser parms(1024);
  parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input file", true, "projections.hdf5" );
  parms.add_parameter( 'o', COMMAND_LINE_INT, 1, "Number of samples per ray", true, "0" );
  parms.add_parameter( 'p', COMMAND_LINE_INT, 1, "Projections per batch", true, "1" );

  parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, 
                       "Output image file name (.real)", true, "reconstruction.real" );
  parms.add_parameter( 'x', COMMAND_LINE_INT,    1, "Matrix size in the x direction", true, "512" );
  parms.add_parameter( 'y', COMMAND_LINE_INT,    1, "Matrix size in the y direction", true, "512" );
  parms.add_parameter( 'z', COMMAND_LINE_INT,    1, "Matrix size in the z direction", true, "1" );

  parms.add_parameter( 'P', COMMAND_LINE_STRING,    1,
                       "Projection space geometry file", true, "ps_geometry.hdf5" );
  parms.add_parameter( 'B', COMMAND_LINE_STRING,    1,
                       "Projection space binning file", true, "binning.hdf5" );
  parms.add_parameter( 'S', COMMAND_LINE_INT, 1, "Use exact SAG correction if present", true, "0" );

  parms.add_parameter( 'X', COMMAND_LINE_STRING,    1, 
                       "Image space spacing in mm - x direction", true, "0.488" );
  parms.add_parameter( 'Y', COMMAND_LINE_STRING,    1, 
                       "Image space spacing in mm - y direction", true, "0.488" );
  parms.add_parameter( 'Z', COMMAND_LINE_STRING,    1,
                       "Image space spacing in mm - z direction", true, "1" );

  parms.add_parameter( 'Q', COMMAND_LINE_FLOAT,    3,
                       "Image space spacing in mm - all three directions", false);


  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ) {
      //logger.info << " Running reconstruction with the following parameters: " << logger.end;
      parms.print_parameter_list();
  }
  else{
    parms.print_parameter_list();
    parms.print_usage();
    return 1;
  }
  
  
  std::string psgeometrydata_filename = (char*)parms.get_parameter('P')->get_string_value();
  std::cout << "ps geometry data file: " << psgeometrydata_filename << std::endl;
  PS_Geometry* ps_g = new PS_Geometry();
  ps_g->loadData(psgeometrydata_filename);
  ps_g->print(std::cout);




  std::string inFile = (char*)parms.get_parameter('d')->get_string_value();
  std::string outFile = (char*)parms.get_parameter('r')->get_string_value();


  PS_Dataset* ps = new PS_Dataset(ps_g);
  ps->loadData(inFile, false, false);
  unsigned int ip_w = ps->getProjections()->get_size(0);
  unsigned int ip_h = ps->getProjections()->get_size(1);
  uintd2 ps_dims_in_pixels =uintd2(ip_w, ip_h);


  float ip_sx = ps_g->getSpacingArray()[0];
  float ip_sy = ps_g->getSpacingArray()[1];

  float ip_x = ip_sx * ps_dims_in_pixels[0];
  float ip_y = ip_sy * ps_dims_in_pixels[1];

  ps_g->setSpacing(ip_sx,ip_sy);
  float SDD = ps_g->getSDD();
  float SAD = ps_g->getSAD();
  // Projection space physical dimensions 40cm x 30cm
  floatd2 ps_dims_in_mm = floatd2(ip_x, ip_y);

  unsigned int is_x_in_pixels = parms.get_parameter('x')->get_int_value();
  unsigned int is_y_in_pixels = parms.get_parameter('y')->get_int_value();
  unsigned int is_z_in_pixels = parms.get_parameter('z')->get_int_value();
  floatd3 is_dims_in_pixels = floatd3(is_x_in_pixels, is_y_in_pixels, is_z_in_pixels); // in pixels

  floatd3 is_dims_in_mm;
  floatd3 is_spacing_in_mm;
  if (parms.get_parameter('Q')->get_is_set()){
	  float tmp_x = parms.get_parameter('Q')->get_float_value(0);
	  float tmp_y = parms.get_parameter('Q')->get_float_value(1);
	  float tmp_z = parms.get_parameter('Q')->get_float_value(2);
      is_dims_in_mm = floatd3(tmp_x, tmp_y, tmp_z);
      is_spacing_in_mm = is_dims_in_mm/is_dims_in_pixels;
  } else {
      float is_spacing_x_in_mm = atof(parms.get_parameter('X')->get_string_value()); // in mm. for full fan
      float is_spacing_y_in_mm = atof(parms.get_parameter('Y')->get_string_value()); // in mm. for full fan
      float is_spacing_z_in_mm = atof(parms.get_parameter('Z')->get_string_value()); // in mm.
      is_spacing_in_mm = floatd3(is_spacing_x_in_mm, is_spacing_y_in_mm, is_spacing_z_in_mm);
      is_dims_in_mm = is_spacing_in_mm * is_dims_in_pixels;
  }


  bool use_sag_correction = parms.get_parameter('S')->get_int_value();

  PS_BinningData* ps_bd = new PS_BinningData();
	ps_bd->generateData(ps_g);
  ps_bd->print(std::cout);

  unsigned int ppb = parms.get_parameter('p')->get_int_value(); // projections per batch
  //assert(ppb==1); // there is an error with ppb other than 1


  unsigned int numSamplesPerRay = parms.get_parameter('o')->get_int_value();


  //floatd3 is_dims_in_mm = is_dims_in_pixels * is_spacing_in_mm;
  float diameter_in_mm = norm( (is_dims_in_mm + is_spacing_in_mm) );
  float lengthOfRay_in_mm = norm(is_dims_in_mm);
  unsigned int numSamplesPerPixel = 3;
  float minSpacing = min(is_spacing_in_mm)/numSamplesPerPixel;
  if (numSamplesPerRay == 0) {
      numSamplesPerRay = ceil( lengthOfRay_in_mm / minSpacing );
  }

  float step_size_in_mm = lengthOfRay_in_mm / numSamplesPerRay;

  if (!use_sag_correction) {
  	ps_g->setSAGx(0,0,0);
  	ps_g->setSAGy(0,0,0);
  }
  size_t numProjs = ps->getProjections()->get_size(2);
  size_t height = ip_h;
  size_t width = ip_w;

  size_t needed_bytes = 2 * prod(is_dims_in_pixels) * sizeof(_real);




  std::vector<unsigned int> is_dims;
  is_dims.push_back(is_dims_in_pixels[0]);
  is_dims.push_back(is_dims_in_pixels[1]);
  is_dims.push_back(is_dims_in_pixels[2]);

  std::cout << "IS dimensions " << is_dims[0] << " " << is_dims[1] << " " << is_dims[2] << std::endl;

  std::vector<unsigned int> ps_dims;
  ps_dims.push_back(ps_dims_in_pixels[0]);
  ps_dims.push_back(ps_dims_in_pixels[1]);
  ps_dims.push_back(numProjs);
  boost::shared_ptr< hoCuNDArray<_real> > projections = boost::static_pointer_cast<hoCuNDArray<_real> >(ps->getProjections());

  //Standard 3d FDK
  // Define encoding matrix
  boost::shared_ptr< hoCudaConebeamProjectionOperator<_real> > 
      E( new hoCudaConebeamProjectionOperator<_real>() );
  E->setup( ps_g, ps_bd, ps_g->getAnglesArray(), ppb,
                 is_spacing_in_mm, ps_dims_in_pixels,
                 numSamplesPerRay,  true);
  E->set_codomain_dimensions(projections->get_dimensions().get());
  // Form right hand side
  E->set_domain_dimensions(&is_dims);

  hoCuNDArray<_real> fdk(&is_dims);
  E->mult_MH(projections.get(),&fdk);


  write_nd_array<_real>( &fdk,"fdk.real" );

  /*4D FDK-MB algorithm starts here. McKinnon GC, RHT Bates,
   *
   *"Towards Imaging the Beating Heart Usefully with a Conventional CT Scanner,"
   *" Biomedical Engineering, IEEE Transactions on , vol.BME-28, no.2, pp.123,127, Feb. 1981
   * doi: 10.1109/TBME.1981.324785
  */



  PS_BinningData* ps_bd4d = new PS_BinningData();
 	std::string binningdata_filename = (char*)parms.get_parameter('B')->get_string_value();
 	std::cout << "binning data file: " << binningdata_filename << std::endl;
 	ps_bd4d->loadData(binningdata_filename);
 	ps_bd4d->print(std::cout);

 	size_t numBins = ps_bd4d->getBinningData().size();
 	is_dims.push_back(numBins);
	boost::shared_ptr< hoCudaConebeamProjectionOperator<_real> >
	      E4D( new hoCudaConebeamProjectionOperator<_real>() );
	E4D->setup( ps_g, ps_bd4d, ps_g->getAnglesArray(), ppb,
								 is_spacing_in_mm, ps_dims_in_pixels,
								 numSamplesPerRay,  true);
	E4D->set_codomain_dimensions(projections->get_dimensions().get());
	// Form right hand side
	E4D->set_domain_dimensions(&is_dims);


	hoCuNDArray<_real> diff_proj(projections->get_dimensions());

	E->mult_M(&fdk,&diff_proj);
	*projections -= diff_proj;

	hoCuNDArray<_real> result(&is_dims);
	E4D->mult_MH(projections.get(),&result);

	result += fdk;




	write_nd_array<_real>( &result, outFile.c_str() );


}
