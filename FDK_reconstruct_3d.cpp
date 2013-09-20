#include "parameterparser.h"
#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCudaConebeamProjectionOperator.h"
#include "hoNDArray_fileio.h"
#include "hoCuNDArray_math.h"
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
	parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input acquisition filename (.hdf5)", true );
	parms.add_parameter( 'B', COMMAND_LINE_STRING, 1, "Binning filename (.hdf5) - 1st phase only", false );
	parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image filename (.real)", true, "reconstruction_FDK.real" );
	parms.add_parameter( 'm', COMMAND_LINE_INT, 3, "Matrix size (3d)", true, "256, 256, 144" );
	parms.add_parameter( 'f', COMMAND_LINE_FLOAT, 3, "FOV in mm (3d)", true, "448, 448, 252" );
	parms.add_parameter( 'F', COMMAND_LINE_INT, 1, "Use filtered backprojection (fbp)", true, "1" );
	parms.add_parameter( 'O', COMMAND_LINE_INT, 1, "Use oversampling in fbp", false );
	parms.add_parameter( 'P', COMMAND_LINE_INT, 1, "Projections per batch", false );

	parms.parse_parameter_list(argc, argv);
	if( parms.all_required_parameters_set() ) {
		parms.print_parameter_list();
	}
	else{
		parms.print_parameter_list();
		parms.print_usage();
		return 1;
	}

	std::string acquisition_filename = (char*)parms.get_parameter('d')->get_string_value();
	std::string image_filename = (char*)parms.get_parameter('r')->get_string_value();

	// Load acquisition data
	//

	boost::shared_ptr<CBCT_acquisition> acquisition( new CBCT_acquisition() );
	acquisition->load(acquisition_filename);

	/*
	std::vector<floatd2> offsets = acquisition->get_geometry()->get_offsets();
	std::fill(offsets.begin(),offsets.end(),floatd2(0,0));
	acquisition->get_geometry()->set_offsets(offsets);*/
	// Generate default binning data
	//
	boost::shared_ptr<CBCT_binning> binning( new CBCT_binning() );

	if (parms.get_parameter('B')->get_is_set()){
			std::string binningdata_filename = (char*)parms.get_parameter('B')->get_string_value();
			std::cout << "binning data file: " << binningdata_filename << std::endl;
			binning->load(binningdata_filename);
	} else binning->set_as_default_3d_bin(acquisition->get_projections()->get_size(2));



	// Configuring...
	//

	uintd2 ps_dims_in_pixels( acquisition->get_projections()->get_size(0),
			acquisition->get_projections()->get_size(1) );

	floatd2 ps_dims_in_mm( acquisition->get_geometry()->get_FOV()[0],
			acquisition->get_geometry()->get_FOV()[1] );

	float SDD = acquisition->get_geometry()->get_SDD();
	float SAD = acquisition->get_geometry()->get_SAD();

	uintd3 is_dims_in_pixels( parms.get_parameter('m')->get_int_value(0),
			parms.get_parameter('m')->get_int_value(1),
			parms.get_parameter('m')->get_int_value(2) );

	floatd3 is_dims_in_mm( parms.get_parameter('f')->get_float_value(0),
			parms.get_parameter('f')->get_float_value(1),
			parms.get_parameter('f')->get_float_value(2) );

	bool use_fbp = parms.get_parameter('F')->get_int_value();

	// Allocate array to hold the result
	//

	std::vector<unsigned int> is_dims;
	is_dims.push_back(is_dims_in_pixels[0]);
	is_dims.push_back(is_dims_in_pixels[1]);
	is_dims.push_back(is_dims_in_pixels[2]);
	is_dims.push_back(1); // one temporal frame for 3d

	hoCuNDArray<float> fdk_3d(&is_dims);

	//
	// Standard 3D FDK reconstruction
	//

	boost::shared_ptr< hoCudaConebeamProjectionOperator > E( new hoCudaConebeamProjectionOperator() );

	E->setup( acquisition, binning, is_dims_in_mm );
	E->use_filtered_backprojections(use_fbp);

	CommandLineParameter *parm = parms.get_parameter('P');
	if( parm && parm->get_is_set() )
		E->set_num_projections_per_batch( parm->get_int_value() );

	parm = parms.get_parameter('O');
	if( parm && parm->get_is_set() )
		E->use_oversampling_in_filtered_backprojection( bool(parm->get_int_value()) );

	{
		GPUTimer timer("Running 3D FDK reconstruction");
		hoCuNDArray<float> projections(*acquisition->get_projections());
		E->mult_MH( &projections, &fdk_3d );
	}

	write_nd_array<float>( &fdk_3d, image_filename.c_str() );
	return 0;
}
