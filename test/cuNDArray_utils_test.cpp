#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_elemwise.h"
#include "NFFTOperator.h"
#include "cuNFFT.h"
#include "complext.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "vector_td_utilities.h"

// Std includes
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "GadgetronTimer.h"
#include "hoArmadillo.h"

using namespace Gadgetron;
using testing::Types;

template <typename T> class cuNDArray_utils_TestReal : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(&dims);
    Array2 = cuNDArray<T>(&dims);
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
};

template <typename T> class cuNDArray_utils_TestCplx : public ::testing::Test {
protected:
  virtual void SetUp() {
    size_t vdims[] = {37, 49, 23, 19}; //Using prime numbers for setup because they are messy
    dims = std::vector<size_t>(vdims,vdims+sizeof(vdims)/sizeof(size_t));
    Array = cuNDArray<T>(&dims);
    Array2 = cuNDArray<T>(&dims);

    // Prep for NUFFT_TEST
    RO = 2500;
    INT = 1100;
    CHA = 8;
    xsize_ = 256;
    ysize_ = 256;
    kernel_width_ = 3;
    oversampling_factor_ = 2.1;


    std::vector<size_t> data_dims = {RO, INT,CHA};
    fake_data=cuNDArray<float_complext>(data_dims); 
    data_dims.pop_back();
    fake_dcw=cuNDArray<float>(data_dims); 
    
    hoNDArray<floatd2> fake_traj_ho(data_dims);
    vector_td<float, 2> init_val;
    init_val[0]= 0.1f;
    init_val[1]= 0.1f;
    fake_traj_ho.fill(init_val);

    fake_traj = cuNDArray<floatd2>(fake_traj_ho);
    fill(&fake_dcw,1.0f);
    fill(&fake_data,complext(1.0f,0.0f));

    
    unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();

    image_dims_.push_back(xsize_);
    image_dims_.push_back(ysize_);

    image_dims_os_ = uint64d2(((static_cast<size_t>(std::ceil(image_dims_[0] * oversampling_factor_)) + warp_size - 1) / warp_size) * warp_size,
                              ((static_cast<size_t>(std::ceil(image_dims_[1] * oversampling_factor_)) + warp_size - 1) / warp_size) * warp_size); // No oversampling is needed in the z-direction for SOS

    recon_dims = {this->image_dims_[0], this->image_dims_[1], CHA};

   
  }
  std::vector<size_t> dims;
  cuNDArray<T> Array;
  cuNDArray<T> Array2;
  
  cuNDArray<float_complext> fake_data;
  cuNDArray<float>   fake_dcw;
  cuNDArray<floatd2> fake_traj;
  size_t RO, INT, Kz, CHA, xsize_, ysize_; 
  float kernel_width_, oversampling_factor_;
  uint64d2 image_dims_os_;
  std::vector<size_t> image_dims_;
  boost::shared_ptr<cuNFFT_plan<float, 2>> nfft_plan_;
  std::vector<size_t> recon_dims;
  
};

typedef Types<float, double> realImplementations;
typedef Types</*std::complex<float>, std::complex<double>,*/ float_complext, double_complext> cplxImplementations;

TYPED_TEST_SUITE(cuNDArray_utils_TestReal, realImplementations);

TYPED_TEST(cuNDArray_utils_TestReal,permuteTest){

  fill(&this->Array,TypeParam(1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  TypeParam tmp(2);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, permute(this->Array,order).at(0));
  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(37));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(1));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(19));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, permute(this->Array,order).at(851));
}

TYPED_TEST(cuNDArray_utils_TestReal,shiftDimTest){

  fill(&this->Array,TypeParam(1));

  TypeParam tmp(2);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, shift_dim(this->Array,0).at(0));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,0).at(37));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,1).at(1));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,-1).at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,2).at(23*37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,3).at(37*19));
  EXPECT_FLOAT_EQ(2, shift_dim(this->Array,4).at(37));
}

TYPED_TEST(cuNDArray_utils_TestReal,sumTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(49*v1,sum(&this->Array,1)->at(idx));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(23*v1,sum(&this->Array,2)->at(idx));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(19*v1,sum(&this->Array,3)->at(idx));
}


TYPED_TEST(cuNDArray_utils_TestReal,meanTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_NEAR(v1,mean(&this->Array), 0.001);

}
TYPED_TEST_SUITE(cuNDArray_utils_TestCplx, cplxImplementations);



TYPED_TEST(cuNDArray_utils_TestCplx,meanTest){
  TypeParam v1 = TypeParam(12.34);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(v1),real(mean(&this->Array)));

}

TYPED_TEST(cuNDArray_utils_TestCplx,permuteTest){
  
  fill(&this->Array,TypeParam(1,1));

  std::vector<size_t> order;
  order.push_back(0); order.push_back(1); order.push_back(2); order.push_back(3);
  
  TypeParam tmp(2,3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, real(permute(this->Array,order).at(0)));
  EXPECT_FLOAT_EQ(1, imag(permute(this->Array,order).at(0)));

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(37)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(37)));

  order.clear();
  order.push_back(1); order.push_back(0); order.push_back(2); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(1)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(1)));

  order.clear();
  order.push_back(3); order.push_back(1); order.push_back(2); order.push_back(0);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(19)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(19)));

  order.clear();
  order.push_back(2); order.push_back(0); order.push_back(1); order.push_back(3);

  EXPECT_FLOAT_EQ(2, real(permute(this->Array,order).at(851)));
  EXPECT_FLOAT_EQ(3, imag(permute(this->Array,order).at(851)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,shiftDimTest){

  fill(&this->Array,TypeParam(1,1));

  TypeParam tmp(2,3);
  CUDA_CALL(cudaMemcpy(&this->Array.get_data_ptr()[37], &tmp, sizeof(TypeParam), cudaMemcpyHostToDevice));

  EXPECT_FLOAT_EQ(1, real(shift_dim(this->Array,0).at(0)));
  EXPECT_FLOAT_EQ(1, imag(shift_dim(this->Array,0).at(0)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,0).at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,0).at(37)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,1).at(1)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,1).at(1)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,-1).at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,-1).at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,2).at(23*37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,2).at(23*37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,3).at(37*19)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,3).at(37*19)));

  EXPECT_FLOAT_EQ(2, real(shift_dim(this->Array,4).at(37)));
  EXPECT_FLOAT_EQ(3, imag(shift_dim(this->Array,4).at(37)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,sumTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(49)*v1),real(sum(&this->Array,1)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(49)*v1),imag(sum(&this->Array,1)->at(idx)));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(23)*v1),real(sum(&this->Array,2)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(23)*v1),imag(sum(&this->Array,2)->at(idx)));

  fill(&this->Array,v1);
  EXPECT_FLOAT_EQ(real(TypeParam(19)*v1),real(sum(&this->Array,3)->at(idx)));
  EXPECT_FLOAT_EQ(imag(TypeParam(19)*v1),imag(sum(&this->Array,3)->at(idx)));
}

TYPED_TEST(cuNDArray_utils_TestCplx,padTest){
  TypeParam v1 = TypeParam(12.34, 56.78);
  unsigned int idx = 0;

  fill(&this->Array,v1);

  vector_td<size_t,4> size = from_std_vector<size_t,4>(this->dims);
  size *= 2;

  auto out = pad<TypeParam,4>(size,this->Array);

  double scale = std::pow(2.0,4);
  EXPECT_EQ(out.get_number_of_elements(),this->Array.get_number_of_elements()*scale);
  EXPECT_NEAR(real(mean(&out))*scale,real(mean(&this->Array)), 0.001);
  EXPECT_NEAR(imag(mean(&out))*scale,imag(mean(&this->Array)), 0.001);
}

TYPED_TEST(cuNDArray_utils_TestCplx,cuNFFT_ATOMIC){
        std::vector<size_t> flat_dims = {this->fake_traj.get_number_of_elements()};
     cuNDArray<floatd2> flat_traj(flat_dims, this->fake_traj.get_data_ptr());
    this->nfft_plan_ = NFFT<cuNDArray, float, 2>::make_plan(from_std_vector<size_t, 2>(this->image_dims_), this->image_dims_os_, this->kernel_width_, ConvolutionType::ATOMIC);
    this->nfft_plan_->preprocess(flat_traj, NFFT_prep_mode::NC2C);

     auto temp = boost::make_shared<cuNDArray<float_complext>>(this->recon_dims);

     this->nfft_plan_->compute(&this->fake_data, *temp, &this->fake_dcw, NFFT_comp_mode::BACKWARDS_NC2C);

}

TYPED_TEST(cuNDArray_utils_TestCplx,cuNFFT_STANDARD){
  
    std::vector<size_t> flat_dims = {this->fake_traj.get_number_of_elements()};
     cuNDArray<floatd2> flat_traj(flat_dims, this->fake_traj.get_data_ptr());
     this->nfft_plan_ = NFFT<cuNDArray, float, 2>::make_plan(from_std_vector<size_t, 2>(this->image_dims_), this->image_dims_os_, this->kernel_width_, ConvolutionType::STANDARD);
     this->nfft_plan_->preprocess(flat_traj, NFFT_prep_mode::NC2C);

     
     auto temp = boost::make_shared<cuNDArray<float_complext>>(this->recon_dims);
 
     this->nfft_plan_->compute(&this->fake_data, *temp, &this->fake_dcw, NFFT_comp_mode::BACKWARDS_NC2C);
 
}


TEST(padTest,largeSize){
// So, this test is mainly here because pad apparently fails for large sized arrays.
	size_t vdims[] = {192,192,50};
	std::vector<size_t> dims(vdims,vdims+sizeof(vdims)/sizeof(size_t));
	size_t vdims2[] = {256,256,256};
	std::vector<size_t> dims2(vdims2,vdims2+sizeof(vdims2)/sizeof(size_t));

	cuNDArray<float_complext> in(&dims);
	fill(&in,float_complext(1));
	cuNDArray<float_complext> out(&dims2);

	pad<float_complext,3>(in,out);

	EXPECT_FLOAT_EQ(nrm2(&in),nrm2(&out));

}
