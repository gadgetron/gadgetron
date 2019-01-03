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
#include "hoCuCgSolver.h"
#include "hoCuEncodingOperatorContainer.h"

#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCuConebeamProjectionOperator.h"
#include "hoCuOFConebeamProjectionOperator.h"

#include <string>
#include <sstream>

using namespace Gadgetron;

boost::shared_ptr< hoNDArray<float> >
perform_registration( boost::shared_ptr< hoNDArray<float> > volume, unsigned int phase, float of_alpha, float of_beta, unsigned int num_multires_levels )
{
	std::vector<size_t> volume_dims_3d = *volume->get_dimensions();
	volume_dims_3d.pop_back();
	std::vector<size_t> volume_dims_3d_3 = volume_dims_3d;
	volume_dims_3d_3.push_back(3);

	size_t num_elements_3d = volume_dims_3d[0]* volume_dims_3d[1]* volume_dims_3d[2];
	size_t num_phases = volume->get_size(3);

	boost::shared_ptr< hoNDArray<float> > host_result_field( new hoNDArray<float> );
	{
		std::vector<size_t> volume_dims_4d = volume_dims_3d;
		volume_dims_4d.push_back(3);
		volume_dims_4d.push_back(volume->get_size(3)-1);
		host_result_field->create( &volume_dims_4d );
	}

	// Upload host data to device
	//

	unsigned int counter = 0;
	hoNDArray<float> host_moving( &volume_dims_3d, volume->get_data_ptr()+phase*num_elements_3d );

	for( unsigned int i=0; i<num_phases; i++ ){

		if( i==phase )
			continue;

		hoNDArray<float> host_fixed( &volume_dims_3d, volume->get_data_ptr()+i*num_elements_3d );

		cuNDArray<float> fixed_image(&host_fixed);
		cuNDArray<float> moving_image(&host_moving);

		boost::shared_ptr< cuLinearResampleOperator<float,3> > R( new cuLinearResampleOperator<float,3>() );

		// Setup solver
		//

		cuCKOpticalFlowSolver<float,3> OFs;
		OFs.set_interpolator( R );
		OFs.set_output_mode( cuCKOpticalFlowSolver<float,3>::OUTPUT_VERBOSE );
		OFs.set_max_num_iterations_per_level( 500 );
		OFs.set_num_multires_levels( num_multires_levels );
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
		std::vector<size_t> order;
		order.push_back(0); order.push_back(1); order.push_back(2);
		order.push_back(4); order.push_back(3);
		cuNDArray<float> tmp(host_result_field.get()); // permute is too slow on the host
		host_result_field = permute(tmp, order).to_host();

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
	//parms.add_parameter( 'T', COMMAND_LINE_FLOAT, 1, "TV regularization weight", true, "0.0" );

	parms.add_parameter( 'P', COMMAND_LINE_INT, 1, "Projections per batch", false );
	parms.add_parameter( 'D', COMMAND_LINE_INT, 1, "Number of downsamples of projection plate", true, "0" );
	parms.add_parameter( 'R', COMMAND_LINE_INT, 1, "Number of downsamples of registration field", true, "1" );

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

	boost::shared_ptr<CBCT_binning> binning_4d( new CBCT_binning() );
	{
		GPUTimer timer("Loading binning data");
		binning_4d->load(binning_filename);
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

	size_t num_phases = binning_4d->get_number_of_bins();

	float of_alpha = parms.get_parameter('A')->get_float_value();
	float of_beta = parms.get_parameter('B')->get_float_value();
	//float tv_weight = parms.get_parameter('T')->get_float_value();

	std::vector<size_t> is_dims_3d;
	is_dims_3d.push_back(is_dims_in_pixels[0]);
	is_dims_3d.push_back(is_dims_in_pixels[1]);
	is_dims_3d.push_back(is_dims_in_pixels[2]);

	unsigned int num_elements_3d = is_dims_in_pixels[0]*is_dims_in_pixels[1]*is_dims_in_pixels[2];

	std::vector<size_t> is_dims_4d = is_dims_3d;
	is_dims_4d.push_back(num_phases);

  // Create mask to zero-out the areas of the images that ar enot fully sampled
  // Do this by a 3D unfiltered backprojection of projections of pure ones.
  // - threshold on a 2% level
  //
    
  hoCuNDArray<float> image_filter( &is_dims_3d );

  {
		GPUTimer timer("Filtering TV reconstruction");
    hoCuNDArray<float> projections( *acquisition->get_projections() );
    boost::shared_ptr<CBCT_binning> binning_3d( new CBCT_binning );
    hoCuConebeamProjectionOperator E;
    fill( &projections, 1.0f );
    binning_3d->set_as_default_3d_bin(projections.get_size(2));
    E.setup( acquisition, binning_3d, is_dims_in_mm );
    E.set_use_filtered_backprojection(false);
    E.mult_MH( &projections, &image_filter ); 
    if( !E.get_use_offset_correction() ){
      std::cout << std::endl << "Error: currently offset correction is assumed to guide the TV filtering. Offset correction unavailable" << std::endl;
      exit(1);
    }

    float threshold = 0.98f * 0.5f; // 0.5 is the intensity after offset correction
    clamp( &image_filter, threshold, threshold, 0.0f, 1.0f );

    *tv_recon *= image_filter;
  }

  // Downsample TV recon if specified
  // - for every downsample of OF image space, dowmsample OF projections as well
  //
  
  unsigned int num_reg_downsamples = parms.get_parameter('R')->get_int_value();

  boost::shared_ptr<hoCuNDArray<float> > projections( new hoCuNDArray<float>( *acquisition->get_projections() ));
  boost::shared_ptr<hoCuNDArray<float> > OF_projections = projections;

  if( num_reg_downsamples > 0 ) {
    GPUTimer timer("Downsampling TV reconstruction (and OF operator projections accordingly)");

    std::vector<size_t> tmp_dims_3d = is_dims_3d;
    for( unsigned int i=0; i<tmp_dims_3d.size(); i++ ){
      for( unsigned int d=0; d<num_reg_downsamples; d++ ){
        if( (tmp_dims_3d[i]%2)==1 )
          throw std::runtime_error("Error: input volume for registration must have even size in all dimensions (at all levels) in order to downsample");
        tmp_dims_3d[i] /= 2;
      }
    }
    
    cuNDArray<float> tmp_image_in(tv_recon.get());
    cuNDArray<float> tmp_proj_in(OF_projections.get());
      
    std::vector<size_t> volume_dims_4d = *tv_recon->get_dimensions();
    std::vector<size_t> proj_dims_3d = *projections->get_dimensions();
    
    for( unsigned int d=0; d<num_reg_downsamples; d++ ){
      
      for( unsigned int i=0; i<3; i++ ) volume_dims_4d[i] /= 2; // do not downsample temporal dimension
      for( unsigned int i=0; i<2; i++ ) proj_dims_3d[i] /= 2;   // do not downsample #projections dimension
      
      cuNDArray<float> tmp_image_out(&volume_dims_4d);      
      cuNDArray<float> tmp_proj_out(&proj_dims_3d);
      
      cuDownsampleOperator<float,3> D_image;
      cuDownsampleOperator<float,2> D_proj;
      
      D_image.mult_M( &tmp_image_in, &tmp_image_out );
      D_proj.mult_M( &tmp_proj_in, &tmp_proj_out );
      
      tmp_image_in = tmp_image_out;
      tmp_proj_in = tmp_proj_out;
    }
    
    tv_recon = tmp_image_in.to_host();
    OF_projections = boost::shared_ptr<hoCuNDArray<float> >( new hoCuNDArray<float>( *tmp_proj_in.to_host() ));
  }
  
	// Allocate 3d array for phase-by-phase reconstruction
	// and 4D array to hold the overall result
	//

	hoCuNDArray<float> image_3d(&is_dims_3d);
	hoNDArray<float> image_4d(&is_dims_4d);

	// Define encoding operator for the reconstruction -- a "plain" CBCT operator
	//

	boost::shared_ptr< hoCuConebeamProjectionOperator > E( new hoCuConebeamProjectionOperator() );
	E->setup( acquisition, is_dims_in_mm );
	E->set_domain_dimensions(&is_dims_3d);
	E->set_codomain_dimensions(acquisition->get_projections()->get_dimensions().get());

	// Define the optical flow regularization operator 
	//
  
	boost::shared_ptr< hoCuOFConebeamProjectionOperator > OF( new hoCuOFConebeamProjectionOperator() );
  OF->setup( acquisition, binning_4d, is_dims_in_mm );
	OF->set_domain_dimensions(&is_dims_3d);
	OF->set_codomain_dimensions(OF_projections->get_dimensions().get());
	//OF->set_weight(1.0f/float(num_phases));

  // Combine the two operators in an operator container
  //
  
  boost::shared_ptr< hoCuEncodingOperatorContainer<float> > opContainer( new hoCuEncodingOperatorContainer<float>() );
  
  opContainer->add_operator(E);
  opContainer->add_operator(OF);

  std::vector<hoCuNDArray<float>*> codoms;
  codoms.push_back( projections.get() );
  codoms.push_back( OF_projections.get() );

  boost::shared_ptr< hoCuNDArray<float> > f = opContainer->create_codomain( codoms );
  codoms.clear(); projections.reset(); OF_projections.reset();

	// Setup solver
	//

  hoCuCgSolver<float> solver;
	//hoCuNCGSolver<float> solver;

	solver.set_encoding_operator( opContainer );
  
	/*if( tv_weight > 0.0f ){
		boost::shared_ptr<hoCuTvOperator<float,3> > tvOp(new hoCuTvOperator<float,3>());
		tvOp->set_weight(tv_weight);
		tvOp->set_limit(float(1e-6));
		solver.add_nonlinear_operator(tvOp);
    }*/

	//solver.set_non_negativity_constraint(true);
	//solver.set_output_mode(hoCuNCGSolver<float>::OUTPUT_VERBOSE);
  solver.set_output_mode(hoCuCgSolver<float>::OUTPUT_VERBOSE);
	solver.set_max_iterations(parms.get_parameter('i')->get_int_value());
	solver.set_tc_tolerance(float(1e-9));

	// Solve
	//

	for( unsigned int phase=0; phase<binning_4d->get_number_of_bins(); phase++ ){

    // Define the binning for the current phase
    //

    std::vector<unsigned int> bin = binning_4d->get_bin(phase);
    boost::shared_ptr<CBCT_binning> binning_phase( new CBCT_binning() );
    binning_phase->set_bin( bin, 0 );
    E->set_binning(binning_phase);

    OF->set_encoding_phase(phase);

    boost::shared_ptr< hoNDArray<float> > displacements =
      perform_registration( tv_recon, phase, of_alpha, of_beta, 3 - num_reg_downsamples );
    
    OF->set_displacement_field(displacements);
    displacements.reset();

		boost::shared_ptr<hoNDArray<float> > result_phase = solver.solve( f.get() );

		// Copy result to 4d array
		//

		hoNDArray<float> host_result_3d( &is_dims_3d, image_4d.get_data_ptr()+phase*num_elements_3d );
		host_result_3d = *result_phase;

    // Apply "cropping" filter
    //

    host_result_3d *= image_filter;


		// Write out every 3d volume
		//

		/*{
			char filename[256];
			sprintf(&(filename[0]), "reconstruction_3d_%i.real", phase);
			write_nd_array<float>( result_phase.get(), (char*)filename );
			}*/
	}

	write_nd_array<float>( &image_4d, result_filename.c_str() );

	return 0;
}
