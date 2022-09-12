/**
	\brief command line tool for using the cpu hoNFFT toolbox

	Handles 2D data with the weights provided

	\param d: the kspace data points
	\param k: the kspace trajectory
	\param w: density compensation
	\param s: oversampling factor for gridding recon
	\param q: width of the kernl on oversampled grid
	\param n: non-oversampled size of image
	\param o: output file
*/

#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include <iostream>
#include "parameterparser.h"
#include "vector_td.h"
#include "hoNFFT.h"
#include "hoCgSolver.h"
#include <boost/make_shared.hpp>
#include "hoLsqrSolver.h"

using namespace std;
using namespace Gadgetron;

int main(int argc, char** argv){
	ParameterParser parms;
	parms.add_parameter('d', COMMAND_LINE_STRING, 1, "Data (.cplx)", true, "d.cplx");
	parms.add_parameter('k', COMMAND_LINE_STRING, 1, "Trajectories (.real", true, "k.real");
	parms.add_parameter('w', COMMAND_LINE_STRING, 1, "Density compensation (.real)", true, "w.real");
	parms.add_parameter('s', COMMAND_LINE_FLOAT, 1, "Oversampling factor (float)", true, "1.5");
	parms.add_parameter('q', COMMAND_LINE_FLOAT, 1, "Kernel width (float)", true, "7");
	parms.add_parameter('n', COMMAND_LINE_INT, 1, "Image size (int)", true, "128");
	parms.add_parameter('o', COMMAND_LINE_STRING, 1, "Output file (.cplx)", true, "out.cplx");

	parms.parse_parameter_list(argc, argv);
	if(parms.all_required_parameters_set()){
		cout << "Running nfft with the following parameter:" << endl;
		parms.print_parameter_list();
	}else{
		cout << "Some required parameters are missing: " << endl;
		parms.print_parameter_list();
		parms.print_usage();
		return 1;
	}

	// Load data from disk
	boost::shared_ptr<hoNDArray<complext<float>>> data = 
		read_nd_array<complext<float>>((char*) parms.get_parameter('d')->get_string_value());
	boost::shared_ptr<hoNDArray<float>> traj = 
		read_nd_array<float>((char*) parms.get_parameter('k')->get_string_value());
	boost::shared_ptr<hoNDArray<float>> weights = 
		read_nd_array<float>((char*) parms.get_parameter('w')->get_string_value());
	float osf = parms.get_parameter('s')->get_float_value();
	float kernelWidth = parms.get_parameter('q')->get_float_value();
	int n = parms.get_parameter('n')->get_int_value();

	// Print data info
	cout << "data n: " << data->get_number_of_elements() << endl;
	cout << "traj n: " << traj->get_number_of_elements() << endl;
	cout << "weights n: " << weights->get_number_of_elements() << endl;
	cout << "osf: " << osf << endl;
	cout << "Kernel width: " << kernelWidth << endl;
	cout << "Matrix size: " << n << endl;

	// Reformat data for hoNFFT
	hoNDArray<typename reald<float,2>::Type> k(weights->get_number_of_elements());
	for(size_t i = 0; i < traj->get_number_of_elements()/2; i++)
		k[i][0] = (*traj)[i];
	for(size_t i = traj->get_number_of_elements()/2; i < traj->get_number_of_elements(); i++)
		k[i-traj->get_number_of_elements()/2][1] = (*traj)[i];
	
	typename uint64d<2>::Type matrixSize; matrixSize[0] = n; matrixSize[1] = n;
	typename uint64d<2>::Type matrixSizeOs; matrixSizeOs[0] = n*osf; matrixSizeOs[1] = n*osf;
	
	hoNFFT_plan<float, 2> plan(matrixSize, osf, kernelWidth);
	hoNDArray<complext<float>> result(matrixSizeOs[0], matrixSizeOs[1]);
	plan.preprocess(k);
	plan.compute((*data), result, (*weights), hoNFFT_plan<float, 2>::NFFT_BACKWARDS_NC2C);
		
	auto output = boost::make_shared<hoNDArray<complext<float>>>(result);
	write_nd_array<complext<float>>(output.get(), (char*) parms.get_parameter('o')->get_string_value());
}

