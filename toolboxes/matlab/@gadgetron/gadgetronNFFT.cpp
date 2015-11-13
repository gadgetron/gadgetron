#include "cuNFFTOperator.h"

//#include "mex.h"
#include "MatlabUtils.h"
#include <boost/make_shared.hpp>
#include "vector_td_operators.h"
#include "cuNDArray_fileio.h"
#include "cuCgSolver.h"
namespace Gadgetron{
template<unsigned int N> static boost::shared_ptr<hoNDArray<float_complext> > gadgetronNFFT_instance(hoNDArray<float_complext> * input_data, hoNDArray<vector_td<float,N> >* trajectory,
		vector_td<uint64_t,N> matrix_size, float W, hoNDArray<float>* dcw = nullptr){

	cuNDArray<float_complext> cuInput(*input_data);
	cuNDArray<vector_td<float,N> > cu_traj(*trajectory);
	auto op = boost::make_shared<cuNFFTOperator<float,N>>();
	op->setup(matrix_size,matrix_size*size_t(2),W);
	op->preprocess(&cu_traj);
	if (dcw){
		auto cu_dcw = boost::make_shared<cuNDArray<float>>(*dcw);
		sqrt_inplace(cu_dcw.get());
		op->set_dcw(cu_dcw);

		cuInput *= *cu_dcw;
	}
	std::vector<size_t> out_dims(&matrix_size[0],&matrix_size[N]);
	out_dims.push_back(cuInput.get_number_of_elements()/cu_traj.get_number_of_elements());
/*
	op->set_domain_dimensions(&out_dims);
	op->set_codomain_dimensions(cuInput.get_dimensions().get());
	cuCgSolver<float_complext> cg;
	cg.set_max_iterations(10);
	cg.set_tc_tolerance(1e-8);
	cg.set_encoding_operator(op);
	auto output = cg.solve(&cuInput);
*/
	cuNDArray<float_complext> output(out_dims);
	op->mult_MH(&cuInput,&output);
	return output.to_host();
}

static mxArray* gadgetronNFFT_internal(mxArray* input_data,mxArray* input_trajectory, mxArray* dimensions, float W, mxArray* dcw=nullptr){

	auto g_input_data = MatlabToHoNDArray<float_complext>(input_data);
	auto g_dimensions = MatlabToHoNDArray<size_t>(dimensions);
	boost::shared_ptr<hoNDArray<float>> g_dcw;
	if (dcw) g_dcw = boost::make_shared<hoNDArray<float>>(MatlabToHoNDArray<float>(dcw));
	auto traj_dims = mxGetDimensions(input_trajectory);
	boost::shared_ptr<hoNDArray<float_complext> > output;
	if (traj_dims[0] == 2){
		auto g_traj = MatlabToHoNDArray<vector_td<float,2 >>(input_trajectory);
		vector_td<uint64_t,2> matrix_size((g_dimensions)[0],(g_dimensions)[1]);
		output = gadgetronNFFT_instance(&g_input_data,&g_traj,matrix_size,W,g_dcw.get());
	} else if (traj_dims[1]){
		auto g_traj = MatlabToHoNDArray<vector_td<float,3 >>(input_trajectory);
		vector_td<uint64_t,3> matrix_size((g_dimensions)[0],(g_dimensions)[1],(g_dimensions)[2]);
		output = gadgetronNFFT_instance(&g_input_data,&g_traj,matrix_size,W,g_dcw.get());
	}

	return hoNDArrayToMatlab(output.get());


}
}

extern"C" void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){


	//mexPrintf("Pie!");

	if (nrhs < 4) return;
	if (nlhs < 1) return;

	mxArray* input_data = prhs[0];
	mxArray* input_trajectory = prhs[1];
	mxArray* dimensions = prhs[2];
	float W = mxGetScalar(prhs[3]);

	mxArray* dcw = nullptr;
	if (nrhs == 5) dcw = prhs[4];

	plhs[0] = Gadgetron::gadgetronNFFT_internal(input_data,input_trajectory,dimensions,W,dcw);


}

