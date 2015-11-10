#include "cuNFFTOperator.h"

#include "mex.h"
#include "MatlabUtils.h"
#include <boost/make_shared.hpp>
using namespace Gadgetron;


template<unsigned int N> static boost::shared_ptr<hoNDArray<float_complext> > gadgetronNFFT_instance(hoNDArray<float_complext> & input_data, hoNDArray<vector_td<float,N> >& trajectory,
    vector_td<uint64_t,N> matrix_size, float W, hoNDArray<float>* dcw = nullptr){

  cuNDArray<float_complext> cuInput(input_data);
  cuNDArray<vector_td<float,N> > cu_traj(trajectory);
  cuNFFTOperator<float,N> op;
  op.setup(matrix_size,matrix_size*2,W);
  op.preprocess(&cu_traj);
  if (dcw){
    auto cu_dcw = boost::make_shared<cuNDArray<float>>(*dcw);
    op.set_dcw(cu_dcw);
  }
  cuNDArray<float> output(std::vector<size_t>(&matrix_size[0],&matrix_size[N]));
  op.mult_MH(&input_data,&output);
  return output.to_host();
}

static mxArray* gadgetronNFFT_internal(mxArray* input_data,mxArray* input_trajectory, mxArray* dimensions, mxArray* dcw=nullptr){

  auto g_input_data = MatlabToHoNDArray<float_complext>(input_data);
  auto g_dimensions = MatlabToHoNDArray<size_t>(dimensions);
  boost::shared_ptr<hoNDArray<float>> g_dcw;
  if (dcw) g_dcw = boost::make_shared<hoNDArray<float>>(MatlabToHoNDArray<float>(dcw));
  const mwSize* traj_dims = mxGetDimensions(input_trajectory);
  boost::shared_ptr<hoNDArray<float> > output;
  if (traj_dims[0] == 2){
    auto g_traj = MatlabToHoNDArray<vector_td<float,2 >>(input_trajectory);
    vector_td<uint64_t,2> matrix_size((*g_dimensions)[0],(*g_dimensions)[1]);
    output = gadgetronNFFT_instance(g_input_data,g_traj,matrix_size,W,g_dcw.get());
  } else if (traj_dims[1]){
    auto g_traj = MatlabToHoNDArray<vector_td<float,3 >>(input_trajectory);
    vector_td<uint64_t,3> matrix_size((*g_dimensions)[0],(*g_dimensions)[1],(*g_dimensions)[2]);
    output = gadgetronNFFT_instance(g_input_data,g_traj,matrix_size,W,g_dcw.get());
  }

  return hoNDArrayToMatlab(output.get());


}


void gadgetronNFFT(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){

	const mxArray* input_data = prhs[0];
	const mxArray* input_trajectory = prhs[1];
	const mxArray* dimensions = prhs[2];
	const mxArray* dcw = nullptr;
	if (nrhs == 4) dcw = prhs[3];

	plhs[0] = gadgetronNFFT_internal(input_data,input_trajectory,dimensions,dcw);

}
