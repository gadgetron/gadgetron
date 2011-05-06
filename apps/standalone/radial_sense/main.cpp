#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "vector_td.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "cgOperatorNonCartesianSense.h"
#include "cuCG.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "parameterparser.h"

#include <iostream>

using namespace std;

int main(int argc, char** argv)
{

  //
  // Parse command line
  //

  ParameterParser parms;
  parms.add_parameter('d', COMMAND_LINE_STRING,	1, "Sample data file name", true);
  parms.add_parameter('w', COMMAND_LINE_STRING,	1, "Weights file name", true);
  parms.add_parameter('t', COMMAND_LINE_STRING,	1, "Trajecotry file name", true);
  parms.add_parameter('m', COMMAND_LINE_INT,    1, "Matrix size", true);
  parms.add_parameter('o', COMMAND_LINE_INT,    1, "Oversampled Matrix size", true);
  parms.add_parameter('f', COMMAND_LINE_FLOAT,  1, "Projections per frame", false);
  parms.add_parameter('i', COMMAND_LINE_INT,    1, "Number of iterations", true, "10");
  parms.add_parameter('k', COMMAND_LINE_FLOAT,  1, "Kernel width", true, "5.5");
  parms.add_parameter('K', COMMAND_LINE_FLOAT,  1, "Kappa", true, "0.1");

  parms.parse_parameter_list(argc, argv);
  if( parms.all_required_parameters_set() ){
    cout << " Running reconstruction with the following parameters: " << endl;
    parms.print_parameter_list();
  }
  else{
    cout << " Some required parameters are missing: " << endl;
    parms.print_parameter_list();
    parms.print_usage();
    exit(1);
  }
  
  GPUTimer *timer = new GPUTimer("Loading data");
  hoNDArray<float_complext::Type> host_data = read_nd_array<float_complext::Type>((char*)parms.get_parameter('d')->get_string_value());
  hoNDArray<float> host_weights = read_nd_array<float>((char*)parms.get_parameter('w')->get_string_value());
  hoNDArray<float> host_traj = read_nd_array<float>((char*)parms.get_parameter('t')->get_string_value());
  delete timer;

  host_data.squeeze();
  host_weights.squeeze();
  host_traj.squeeze();

  if( !(host_data.get_number_of_dimensions() == 1 || host_data.get_number_of_dimensions() == 2) || 
      !(host_weights.get_number_of_dimensions() == 1 || host_weights.get_number_of_dimensions() == 2 ) ||
      host_traj.get_number_of_dimensions() != 2 ){
    
    printf("\nInput data is not two-dimensional. Quitting!\n");
    exit(1);
  }
  
  if( host_data.get_size(0) != host_weights.get_size(0) ){
    printf("\nInput data dimensions mismatch between sample data and weights. Quitting!\n");
    exit(1);
  }
  
  if( host_data.get_size(0) != host_traj.get_size(1) || host_traj.get_size(0) != 2 ){
    printf("\nInput trajectory data mismatch. Quitting!\n");
    exit(1);
  }

  // Matrix sizes
  uintd2 matrix_size = uintd2(parms.get_parameter('m')->get_int_value(), parms.get_parameter('m')->get_int_value());
  uintd2 matrix_size_os = uintd2(parms.get_parameter('o')->get_int_value(), parms.get_parameter('o')->get_int_value());
  float kernel_width = parms.get_parameter('k')->get_float_value();
  float kappa = parms.get_parameter('K')->get_float_value();
  unsigned int num_iterations = parms.get_parameter('i')->get_int_value();
  
  // Get trajectory dimensions 'float2' style
  timer = new GPUTimer("uploading to device");
  vector<unsigned int> traj_dims = host_traj.get_dimensions(); traj_dims.erase(traj_dims.begin());
  cuNDArray<float> _traj(host_traj); 
  cuNDArray<floatd2::Type> traj; traj.create(traj_dims, (floatd2::Type*)_traj.get_data_ptr());
  
  cuNDArray<float_complext::Type> *data = new cuNDArray<float_complext::Type>(host_data);
  cuNDArray<float> weights(host_weights);
  delete timer;

  unsigned int num_batches = (data->get_number_of_dimensions() == 2) ? data->get_size(1) : 1;

  //
  // Estimate CSM from initial gridding
  //

  timer = new GPUTimer("Reconstructing for CSM");
  vector<unsigned int> image_dims = cuNDA_toVec<2>(matrix_size); 
  if( num_batches > 1 ) image_dims.push_back(num_batches);
  cuNDArray<float_complext::Type> *image = new cuNDArray<float_complext::Type>(); 
  image->create(image_dims);
  uintd2 fixed_dims(0,0);

  NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, fixed_dims, kernel_width );
  plan.preprocess( &traj, NFFT_plan<float,2>::NFFT_PREP_BACKWARDS );
  plan.compute( data, image, &weights, NFFT_plan<float,2>::NFFT_BACKWARDS );
  delete timer;
  
  timer = new GPUTimer("Estimate CSM");
  auto_ptr< cuNDArray<float_complext::Type> > csm = estimate_b1_map<float,2>( image );
  delete timer; 

  plan.wipe(NFFT_plan<float,2>::NFFT_WIPE_ALL);
  delete data;

  // Form rhs
  std::auto_ptr< cuNDArray<float_complext::Type> > rhs = cuNDA_sum<float_complext::Type>( image, 2 );
  delete image;
  
  cgOperatorNonCartesianSense<float,2> E;
  
  E.setup( matrix_size, matrix_size_os, kernel_width );

  if( E.set_csm(csm.get()) < 0) {
    cout << "Failed to set csm on encoding matrix" << endl;
  }
  
  if (E.set_trajectory(&traj) < 0) {
    cout << "Failed to set trajectory on encoding matrix" << endl;
  }

  if (E.set_weights(&weights) < 0) {
    cout << "Failed to set weights on encoding matrix" << endl;
  }

  //  cuCGPrecondWeight<float2> Dm; brug kappa
  //Dm.set_weights(&D_dev);

  cuCG<float, float_complext::Type> cg;
  cg.add_matrix_operator(&E, 1.0f);
  //cg.set_preconditioner(&Dm);
  cg.set_iterations(num_iterations);
  cg.set_limit(1e-5);
  cg.set_output_mode(cuCG<float, float_complext::Type>::OUTPUT_VERBOSE);

  auto_ptr< cuNDArray<float_complext::Type> > cgresult;
  {
    GPUTimer timer("GPU Conjugate Gradient solve");
    cgresult = cg.solve(rhs.get());
  }

  hoNDArray<float_complext::Type> rho_out = cgresult->to_host();
  write_nd_array<float_complext::Type>(rho_out,"rho_out.cplx");

  hoNDArray<float> host_norm = cuNDA_norm<float,2>(cgresult.get())->to_host();
  write_nd_array<float>( host_norm, "rho_out.real" );

  cout << "Reconstruction done" << endl;
 
  return 0;
}
