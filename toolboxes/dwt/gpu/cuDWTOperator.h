#pragma once

#include "linearOperator.h"
#include "cuNDDWT.h"
#include "vector_td.h"
#include <random>
#include <numeric>
namespace Gadgetron {

template<class T, unsigned int D> class cuDWTOperator : public linearOperator<cuNDArray<T> > {
public:
	cuDWTOperator() : linearOperator<cuNDArray<T> >::linearOperator() {
		run_dimensions = std::vector<size_t>(D,0);
		std::iota(run_dimensions.begin(),run_dimensions.end(),0);
		daubechies4 = vector_td<typename realType<T>::Type, 4>(0.6830127f, 1.1830127f, 0.3169873f, -0.1830127f);
	}

	void set_levels(unsigned int levels_) { levels=levels_;}
	void use_random(bool random){use_random_ = random; }

	virtual ~cuDWTOperator(){};

	virtual void mult_MH_M(cuNDArray<T> * in, cuNDArray<T> * out, bool accumulate = false){
		if (accumulate) *out += *in;
		else *out = *in;
	}
	virtual void mult_M(cuNDArray<T> * in, cuNDArray<T> * out, bool accumulate = false){
		unsigned int loc_levels = levels;
		auto img_dim = in->get_dimensions();
		if (levels == 0 ) loc_levels = calc_levels(img_dim);

		if (use_random_){
			std::default_random_engine generator;
			std::uniform_int_distribution<int> dist(0,3);
			shift_ = dist(generator);
		}
		cuNDArray<T> * tmp_in = in;
		cuNDArray<T> * tmp_out = out;
		if (accumulate ) tmp_out = new cuNDArray<T>(out->get_dimensions());
		if (run_dimensions.size() > 1) tmp_in = new cuNDArray<T>(*in);


		for (auto i = 0; i < loc_levels; i++){
			for (auto dim : run_dimensions){
				cuNDArray<T> small_in(img_dim, tmp_in->get_data_ptr());
				cuNDArray<T> small_out(img_dim, tmp_out->get_data_ptr());
				DWT1<T,D,4>(&small_in,&small_out,daubechies4,dim,shift_);
				std::swap(tmp_in,tmp_out);
			}
			//Resize for next level
			for (int n = 0; n < D; n++)
				img_dim[n] /= 2;
		}

		if (out != tmp_in && !accumulate)*out = *tmp_in;
		if (accumulate)	*out += *tmp_in;
		if (tmp_in != in && tmp_in != out) delete tmp_in;
		if (tmp_out != in && tmp_out != out) delete tmp_out;



	}

	virtual void mult_MH(cuNDArray<T> * in, cuNDArray<T> * out, bool accumulate = false){
		//iint shift = 0;


		cuNDArray<T> * tmp_in = in;
		cuNDArray<T> * tmp_out = out;
		if (accumulate ) tmp_out = new cuNDArray<T>(out->get_dimensions());
		if (run_dimensions.size() > 1) tmp_in = new cuNDArray<T>(*in);

		auto img_dim = in->get_dimensions();
		auto loc_levels = levels;
		if (levels == 0 ) loc_levels = calc_levels(img_dim);
		//Get smallest dimension;
		for (auto n = 0; n < D; n++)
			img_dim[n] /= std::pow(2,loc_levels-1);
		for (auto i = loc_levels; i > 0; i--){
			for (auto dim = run_dimensions.rbegin(); dim != run_dimensions.rend(); ++dim){
				cuNDArray<T> small_in(img_dim, tmp_in->get_data_ptr());
				cuNDArray<T> small_out(img_dim, tmp_out->get_data_ptr());
				IDWT1<T,D,4>(&small_in,&small_out,daubechies4,*dim,shift_);
				std::swap(tmp_in,tmp_out);
			}
			for (auto n = 0; n < D; n++)
				img_dim[n] *= 2;
		}
		if (out != tmp_in && !accumulate)*out = *tmp_in;
		if (accumulate)	*out += *tmp_in;
		if (tmp_in != in && tmp_in != out) delete tmp_in;
		if (tmp_out != in && tmp_out != out) delete tmp_out;

	}

	void set_shift(int shift){ shift_ = shift;}

	vector_td<typename realType<T>::Type, 4> daubechies4;
	/*
	static auto daubechies4 = vector_td<typename realType<T>::Type ,4>{0.6830127f,1.1830127f,0.3169873f,-0.1830127f};
	static auto haahr = vector_td<typename realType<T>::Type,2>{1.0f,1.0f};
	static auto daubechies6= vector_td<typename realType<T>::Type,6>{0.47046721f,1.14111692f,0.650365f,-0.19093442f, -0.12083221f,0.0498175f};
	*/


private:

	unsigned int calc_levels(std::vector<size_t>& dims){
		unsigned int min_dim = std::numeric_limits<unsigned int>::max();
		for (auto dim : run_dimensions){
			min_dim = std::min(min_dim,max_divisions(dims[dim]));
		}
		return min_dim;

	}

	static const unsigned int max_divisions(unsigned int num){
		unsigned int count = 0;
		while (num%2==0 && num > 4) {
			count++;
			num /= 2;
		}
		return count;

	}

	std::vector<size_t> run_dimensions;
	unsigned int levels=0;
	bool use_random_ = false;
	int shift_=0;


};
}
