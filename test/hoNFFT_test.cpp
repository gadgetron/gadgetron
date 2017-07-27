#include "gtest/gtest.h"
#include "complext.h"
#include "hoNDArray_reductions.h"
#include "hoNFFT.h"
#include "vector_td_utilities.h"

using namespace Gadgetron;
using testing::Types;

template<typename Real>
class hoNFFT_2D_NC2C_BACKWARDS : public ::testing::Test{
	protected:
		virtual void SetUp(){
			data.create(10); data.fill(0);
			trajectory.create(10); trajectory.fill(vector_td<Real, 2>(0,0));
			weights.create(10); weights.fill(0.0);
			output.create(30, 30);

			for(int i = 1; i <= 10; i++){
				data[i] = i;
				trajectory[i][0] = Real(i)/Real(10.0)*Real(0.4);
				trajectory[i][1] = Real(i)/Real(10.0)*Real(0.4);
				weights[i] = Real(i)/Real(10.0);
			}
			std::vector<size_t> dims = {20, 20};
			hoNFFT_plan<Real, 2> plan(from_std_vector<size_t, 2>(dims),
				Real(1.5), Real(3));
			plan.preprocess(trajectory);
			plan.compute(data, output, weights, 
				hoNFFT_plan<Real, 2>::NFFT_BACKWARDS_NC2C);
		}
		hoNDArray<complext<Real>> data;
		hoNDArray<vector_td<Real, 2>> trajectory;
		hoNDArray<Real> weights;
		hoNDArray<complext<Real>> output;
};

typedef Types<float, double> realImplementations;

TYPED_TEST_CASE(hoNFFT_2D_NC2C_BACKWARDS, realImplementations);

TYPED_TEST(hoNFFT_2D_NC2C_BACKWARDS, randomTestOne){
	//EXPECT_NEAR(nrm2(&this->output), 1.1503, 1.1503*1e-3);
}
