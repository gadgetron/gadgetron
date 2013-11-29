#include "hoNDArray_fileio.h"
#include "parameterparser.h"
#include "cuNDArray.h"
#include "hoCuNDArray.h"
#include "hoNDArray_utils.h"
#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "hoCuTvOperator.h"
#include "cuCKOpticalFlowSolver.h"
#include "cuLinearResampleOperator.h"
#include "cuDownsampleOperator.h"
#include "GPUTimer.h"
#include "hoCuNCGSolver.h"

#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCuOFProjectionOperator.h"

#include <string>
#include <sstream>

using namespace Gadgetron;

boost::shared_ptr< hoNDArray<float> >
perform_registration( boost::shared_ptr< hoNDArray<float> > volume, unsigned int phase, float of_alpha, float of_beta, bool downsample )
{
	std::vector<unsigned int> volume_dims_3d = *volume->get_dimensions();
	volume_dims_3d.pop_back();

	if( downsample ){
		for( unsigned int i=0; i<volume_dims_3d.size(); i++ ){
			if( (volume_dims_3d[i]%2)==1 )
				throw std::runtime_error("Error: input volume for registration must have even size in all dimensions in order to downsample");
			volume_dims_3d[i] /= 2;
		}
	}

	boost::shared_ptr< hoNDArray<float> > host_volume;

	if( downsample ){
		std::vector<unsigned int> volume_dims_4d = volume_dims_3d;
		volume_dims_4d.push_back(volume->get_dimensions()->at(3));
		cuDownsampleOperator<float,3> D;
		cuNDArray<float> tmp_in(volume.get());
		cuNDArray<float> tmp_out(&volume_dims_4d);
		D.mult_M( &tmp_in, &tmp_out );
		host_volume = tmp_out.to_host();
	}
	else
		host_volume = volume;

	std::vector<unsigned int> volume_dims_3d_3 = volume_dims_3d;
	volume_dims_3d_3.push_back(3);

	unsigned int num_elements_3d = volume_dims_3d[0]* volume_dims_3d[1]* volume_dims_3d[2];
	unsigned int num_phases = host_volume->get_size(3);

	boost::shared_ptr< hoNDArray<float> > host_result_field( new hoNDArray<float> );
	{
		std::vector<unsigned int> volume_dims_4d = volume_dims_3d;
		volume_dims_4d.push_back(3);
		volume_dims_4d.push_back(host_volume->get_size(3)-1);
		host_result_field->create( &volume_dims_4d );
	}

	// Upload host data to device
	//

	unsigned int counter = 0;
	hoNDArray<float> host_moving( &volume_dims_3d, host_volume->get_data_ptr()+phase*num_elements_3d );

	for( unsigned int i=0; i<num_phases; i++ ){

		if( i==phase )
			continue;

		hoNDArray<float> host_fixed( &volume_dims_3d, host_volume->get_data_ptr()+i*num_elements_3d );

		cuNDArray<float> fixed_image(&host_fixed);
		cuNDArray<float> moving_image(&host_moving);

		boost::shared_ptr< cuLinearResampleOperator<float,3> > R( new cuLinearResampleOperator<float,3>() );

		// Setup solver
		//

		cuCKOpticalFlowSolver<float,3> OFs;
		OFs.set_interpolator( R );
		OFs.set_output_mode( cuCKOpticalFlowSolver<float,3>::OUTPUT_VERBOSE );
		OFs.set_max_num_iterations_per_level( 500 );
		OFs.set_num_multires_levels( 3 - ((downsample) ? 1 : 0) );
		OFs.set_alpha(of_alpha);
		OFs.set_beta(of_beta);
		OFs.set_limit(0.01f);

		// Run registration
		//

		boost::shared_ptr< cuNDArray<float> > result = OFs.solve( &fixed_image, &moving_image );

		cuNDArray<float> dev_sub;
		dev_sub.create( &volume_dims_3d_3, result->get_data_ptr() );

		hoNDArray<float> host_sub( &volume_dims_3d_3, host_result_field->get_data_ptr()+counter*num_elements_3d*3 );
		host_sub = *dev_sub.to_host();

		counter++;
	}

	/*
 {
  std::cout << std::endl << "Writing out registration results for phase " << phase << "." << std::endl;
  char filename[256];
  sprintf(&(filename[0]), "def_moving_%i.real", phase);
  write_nd_array<float>(&host_result_image, (char*)filename);
 }
	 */

	{
		// Permute the displacement field (temporal dimension before 'vector' dimension)
		std::vector<unsigned int> order;
		order.push_back(0); order.push_back(1); order.push_back(2);
		order.push_back(4); order.push_back(3);
		cuNDArray<float> tmp(host_result_field.get()); // permute is too slow on the host
		host_result_field = permute(&tmp, &order)->to_host();

		/*char filename[256];
  sprintf(&(filename[0]), "displacement_field_%i.real", phase);
  write_nd_array<float>(host_result_field.get(), (char*)filename);*/
	}

	return host_result_field;
}

int main(int argc, char** argv) 
{
	// Parse command line
	ParameterParser parms(1024);
	parms.add_parameter( 'd', COMMAND_LINE_STRING, 1, "Input file - projections", true, "projections.hdf5" );
	parms.add_parameter( 't', COMMAND_LINE_STRING, 1, "Input file - TV reconstruction", true, "reconstruction_TV.real" );
	parms.add_parameter( 'b', COMMAND_LINE_STRING, 1, "binning file", true, "binning.hdf5" );
	parms.add_parameter( 'r', COMMAND_LINE_STRING, 1, "Output image file name (.real)", true, "reconstruction_REG.real" );

	parms.add_parameter( 'm', COMMAND_LINE_INT, 3, "Matrix size (3d)", true, "256, 256, 144" );
	parms.add_parameter( 'f', COMMAND_LINE_FLOAT, 3, "FOV in mm (3d)", true, "448, 448, 252" );

	parms.add_parameter( 'i', COMMAND_LINE_INT, 1, "Max number of iterations", true, "10" );

	parms.add_parameter( 'A', COMMAND_LINE_FLOAT, 1, "Alpha for OF term", true, "0.05" );
	parms.add_parameter( 'B', COMMAND_LINE_FLOAT, 1, "Beta for OF term", true, "1.0" );
	parms.add_parameter( 'T', COMMAND_LINE_FLOAT, 1, "TV regularization weight", true, "0.0" );

	parms.add_parameter( 'P', COMMAND_LINE_INT, 1, "Projections per batch", false );
	parms.add_parameter( 'D', COMMAND_LINE_INT, 1, "Number of downsamples of projection plate", true, "0" );
	parms.add_parameter( 'R', COMMAND_LINE_INT, 1, "Use downsampling of registration field (bool)", true, "0" );

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
	std::string tv_recon_filename = (char*)parms.get_parameter('t')->get_string_value();
	std::string binning_filename = (char*)parms.get_parameter('b')->get_string_value();
	std::string result_filename = (char*)parms.get_parameter('r')->get_string_value();

	// Load acquisition data
	//

	boost::shared_ptr<CBCT_acquisition> acquisition( new CBCT_acquisition() );
	{
		GPUTimer timer("Loading projections");
		acquisition->load(acquisition_filename);
	}


	// Downsample projections if requested
	//

	{
		GPUTimer timer("Downsampling projections");
		unsigned int num_downsamples = parms.get_parameter('D')->get_int_value();    
		acquisition->downsample(num_downsamples);
	}

	// Load binning data
	//

	boost::shared_ptr<CBCT_binning> binning( new CBCT_binning() );
	{
		GPUTimer timer("Loading binning data");
		binning->load(binning_filename);
	}

	// Load intermediate reconstruction
	// 

	boost::shared_ptr< hoNDArray<float> > tv_recon;
	{
		GPUTimer timer("Loading intermediate reconstruction");
		tv_recon = read_nd_array<float>(tv_recon_filename.c_str());
	}

	if( tv_recon->get_number_of_dimensions() != 4 ){
		printf("\nInput volume for registration must be four-dimensional");
		exit(1);
	}

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

	size_t num_phases = binning->get_number_of_bins();

	float of_alpha = parms.get_parameter('A')->get_float_value();
	float of_beta = parms.get_parameter('B')->get_float_value();
	float tv_weight = parms.get_parameter('T')->get_float_value();

	// Allocate array to hold the result
	//

	std::vector<unsigned int> is_dims_3d;
	is_dims_3d.push_back(is_dims_in_pixels[0]);
	is_dims_3d.push_back(is_dims_in_pixels[1]);
	is_dims_3d.push_back(is_dims_in_pixels[2]);

	unsigned int num_elements_3d = is_dims_in_pixels[0]*is_dims_in_pixels[1]*is_dims_in_pixels[2];

	std::vector<unsigned int> is_dims_4d = is_dims_3d;
	is_dims_4d.push_back(num_phases);

	// Allocate 3d array for phase-by-phase reconstruction
	// and 4D array to hold the overall result
	//

	hoCuNDArray<float> image_3d(&is_dims_3d);
	hoNDArray<float> image_4d(&is_dims_4d);

	bool use_reg_downsampling = bool(parms.get_parameter('R')->get_int_value());

	// Define encoding operator for the reconstruction
	//

	boost::shared_ptr< hoCuOFProjectionOperator > E( new hoCuOFProjectionOperator() );
	E->setup( acquisition, binning, is_dims_in_mm );
	E->set_domain_dimensions(&is_dims_3d);
	E->set_codomain_dimensions(acquisition->get_projections()->get_dimensions().get());

	// Setup solver
	//

	hoCuNCGSolver<float> solver;
	solver.set_encoding_operator(E);

	if( tv_weight > 0.0f ){
		boost::shared_ptr<hoCuTvOperator<float,3> > tvOp(new hoCuTvOperator<float,3>());
		tvOp->set_weight(tv_weight);
		tvOp->set_limit(float(1e-6));
		solver.add_nonlinear_operator(tvOp);
	}

	solver.set_non_negativity_constraint(true);
	solver.set_output_mode(hoCuNCGSolver<float>::OUTPUT_VERBOSE);
	solver.set_max_iterations(parms.get_parameter('i')->get_int_value());
	solver.set_tc_tolerance(float(1e-9));

	// Solve
	//

	hoCuNDArray<float> projections( *acquisition->get_projections() );

	for( unsigned int phase=0; phase<binning->get_number_of_bins(); phase++ ){
		{
			boost::shared_ptr< hoNDArray<float> > displacements =
					perform_registration( tv_recon, phase, of_alpha, of_beta, use_reg_downsampling );

			E->set_encoding_phase(phase);
			E->set_displacement_field(displacements);
		}

		boost::shared_ptr<hoNDArray<float> > result_phase = solver.solve(&projections);

		// Write out every 3d volume
		//

		{
			char filename[256];
			sprintf(&(filename[0]), "reconstruction_3d_%i.real", phase);
			write_nd_array<float>( result_phase.get(), (char*)filename );
		}


		// Copy result to 4d array
		//

		hoNDArray<float> host_result_3d( &is_dims_3d, image_4d.get_data_ptr()+phase*num_elements_3d );
		host_result_3d = *result_phase;
	}

	write_nd_array<float>( &image_4d, result_filename.c_str() );

	return 0;
}
