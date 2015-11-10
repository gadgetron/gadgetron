#include "cuNFFT.h"

#include "mex.h"
#include "MatlabUtils.h"
using namespace Gadgetron;
void gadgetronNFFT(mxArray* input_data,mxArray* input_trajectory, mxArray* dimensions, mxArray* output_array, mxArray* dcw=nullptr){


  auto g_input_data = MatlabToHoNDArray<float_complext>(input_data);
  auto g_dimensions = MatlabToHoNDArray<size_t>(dimensions);
  boost::shared_ptr<hoNDArray<float>> g_dcw;
  if (dcw) g_dcw = boost::make_shared<hoNDArray<float>>(MatlabToHoNDArray<float>(dcw));
  const mwSize* traj_dims = mxGetDimensions(input_trajectory);
  boost::shared_ptr<hoNDArray<float> > output;
  if (traj_dims[0] == 2){
    auto g_traj = MatlabToHoNDArray<vector_td<float,2 >>(input_trajectory);
    vector_td<uint64_t,2> matrix_size(dimensions[0],dimensions[1]);
    output = gadgetronNFFT_instance(g_input_data,g_traj,matrix_size,W,g_dcw.get());
  } else if (traj_dims[1]){
    auto g_traj = MatlabToHoNDArray<vector_td<float,3 >>(input_trajectory);
    vector_td<uint64_t,3> matrix_size(dimensions[0],dimensions[1],dimensions[2]);
    output = gadgetronNFFT_instance(g_input_data,g_traj,matrix_size,W,g_dcw.get());
  }


}

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
  op.mult_MH(&input_data,&output_data);
  return output.to_host();
}
