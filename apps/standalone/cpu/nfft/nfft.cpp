// For easy testing of nfft

#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include <iostream>
#include "parameterparser.h"
#include "vector_td.h"
#include "hoNFFT.h"
#include "hoNFFTOperator.h"
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
	parms.add_parameter('q', COMMAND_LINE_FLOAT, 1, "Kernel width (flaot)", true, "7");
	parms.add_parameter('n', COMMAND_LINE_INT, 1, "Image size (int)", true, "128");
	parms.add_parameter('o', COMMAND_LINE_STRING, 1, "Output file (.cplx)", true, "out.cplx");
	parms.add_parameter('m', COMMAND_LINE_INT, 1, "Mode of the nfft", "0");

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
	int mode = parms.get_parameter('m')->get_int_value();

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

	cout << "Number of data dimensions:" << endl;
	auto dims = *(data->get_dimensions());
	for(auto it : dims) cout << it << ", ";
	cout << endl;

	auto E = boost::make_shared<hoNFFTOperator<float,2>>();
	E->setup(matrixSize, osf, kernelWidth);

	hoLsqrSolver<complext<float>> lsqrSolver;
	lsqrSolver.set_tc_tolerance(1e-8);
	lsqrSolver.set_max_iterations(20);
	lsqrSolver.set_output_mode(hoLsqrSolver<complext<float>>::OUTPUT_VERBOSE);

	hoNDArray<complext<float>> fin(matrixSizeOs[0], matrixSizeOs[0]);
	E->set_domain_dimensions(fin.get_dimensions().get());
	E->set_codomain_dimensions(fin.get_dimensions().get());

	lsqrSolver.set_encoding_operator(E);
	E->preprocess(k);

	hoNFFT_plan<float, 2> p(matrixSize, osf, kernelWidth);
	hoNDArray<complext<float>> r(matrixSizeOs[0], matrixSizeOs[0]);
	p.preprocess(k);
	p.compute(*data, r, *weights, hoNFFT_plan<float, 2>::NFFT_BACKWARDS_NC2C);
	cout << "NORM" << endl;
	cout << Gadgetron::nrm2(&r) << endl;

	if(mode == 0){
		//auto res = lsqrSolver.solve(&r);
		lsqrSolver.solve(&fin, &r);
		write_nd_array<complext<float>>(&fin, (char*) parms.get_parameter('o')->get_string_value());
	}
	
	if(mode == 1){
		hoNFFT_plan<float, 2> plan(matrixSize, osf, kernelWidth);
		hoNDArray<complext<float>> result(matrixSizeOs[0], matrixSizeOs[1]);
		plan.preprocess(k);
		plan.compute((*data), result, (*weights), hoNFFT_plan<float, 2>::NFFT_BACKWARDS_NC2C);
		
		auto output = boost::make_shared<hoNDArray<complext<float>>>(result);
		write_nd_array<complext<float>>(output.get(), (char*) parms.get_parameter('o')->get_string_value());
	}
	
	if(mode == 2){
		hoNFFT_plan<float, 2> plan(matrixSize, osf, kernelWidth);
		hoNDArray<complext<float>> result(weights->get_number_of_elements());
		plan.preprocess(k);
		plan.compute((*data), result, (*weights), hoNFFT_plan<float, 2>::NFFT_FORWARDS_C2NC);
		
		auto output = boost::make_shared<hoNDArray<complext<float>>>(result);
		write_nd_array<complext<float>>(output.get(), (char*) parms.get_parameter('o')->get_string_value());
	}

	if(mode == 3){
		hoNFFT_plan<float, 2> plan(matrixSize, osf, kernelWidth);
		hoNDArray<complext<float>> tmp(192,192);
		hoNDArray<complext<float>> result(weights->get_number_of_elements());

		plan.preprocess(k);
		plan.compute((*data),tmp,(*weights),hoNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);
		plan.compute(tmp,result,(*weights),hoNFFT_plan<float,2>::NFFT_FORWARDS_C2NC);
		plan.compute(result,tmp,(*weights),hoNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);
		auto output = boost::make_shared<hoNDArray<complext<float>>>(tmp);
		write_nd_array<complext<float>>(output.get(), (char*) parms.get_parameter('o')->get_string_value());
		
	}
}

