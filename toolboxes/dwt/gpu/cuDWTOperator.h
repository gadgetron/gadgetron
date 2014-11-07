#pragma once

#include "linearOperator.h"
#include "cuNDDWT.h"
#include "vector_td.h"
#include <random>
namespace Gadgetron {

template<class T, unsigned int D> class cuDWTOperator : public linearOperator<cuNDArray<T> > {
public:
	cuDWTOperator() : linearOperator<cuNDArray<T> >::linearOperator() {
		run_dimensions = std::vector<size_t>(D,0);
		std::iota(run_dimensions.begin(),run_dimensions.end(),0);
	}

	void set_levels(unsigned int levels_) { levels=levels_;}
	void use_random(bool random){use_random_ = random; }

	virtual ~cuDWTOperator(){};

	virtual void mult_MH_M(cuNDArray<T> * in, cuNDArray<T> * out, bool accumulate = false){
		if (accumulate) *out += *in;
		else *out = *in;
	}
	virtual void mult_M(cuNDArray<T> * in, cuNDArray<T> * out, bool accumulate = false){
		int shift = 0;
		if (use_random_){
			std::default_random_engine generator;
			std::uniform_int_distribution<int> dist(0,3);
			shift = dist(generator);
		}
		cuNDArray<T> * tmp_in = in;
		cuNDArray<T> * tmp_out = out;
		if (accumulate ) tmp_out = new cuNDArray<T>(out->get_dimensions());
		if (run_dimensions.size() > 1) tmp_in = new cuNDArray<T>(in);
		auto img_dim = *in->get_dimensions();

		for (auto i = 0; i < levels; i++){
			for (auto dim : run_dimensions){
				cuNDArray<T> small_in(img_dim, tmp_in->get_data_ptr());
				cuNDArray<T> small_out(img_dim, tmp_out->get_data_ptr());
				std::cout << "Dimension " << dim << std::endl;
				DWT1<T,D,4>(&small_in,&small_out,daubechies4,dim,shift);
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
		int shift = 0;
		if (use_random_){
			std::default_random_engine generator;
			std::uniform_int_distribution<int> dist(0,3);
			shift = dist(generator);
		}

		cuNDArray<T> * tmp_in = in;
		cuNDArray<T> * tmp_out = out;
		if (accumulate ) tmp_out = new cuNDArray<T>(out->get_dimensions());
		if (run_dimensions.size() > 1) tmp_in = new cuNDArray<T>(in);

		auto img_dim = *in->get_dimensions();
		//Get smallest dimension;
		for (auto n = 0; n < D; n++)
			img_dim[n] /= std::pow(2,levels-1);

		for (auto i = levels; i > 0; i--){
			for (auto dim = run_dimensions.rbegin(); dim != run_dimensions.rend(); ++dim){
				cuNDArray<T> small_in(img_dim, tmp_in->get_data_ptr());
				cuNDArray<T> small_out(img_dim, tmp_out->get_data_ptr());
				IDWT1<T,D,4>(&small_in,&small_out,daubechies4,*dim,shift);
				std::swap(tmp_in,tmp_out);
			}
			for (auto n = 0; n < D; n++)
				img_dim[n] *= 2;
			std::cout << "Level " << i << std::endl;
		}
		if (out != tmp_in && !accumulate)*out = *tmp_in;
		if (accumulate)	*out += *tmp_in;
		if (tmp_in != in && tmp_in != out) delete tmp_in;
		if (tmp_out != in && tmp_out != out) delete tmp_out;

	}

	constexpr static auto daubechies4 = vector_td<typename realType<T>::Type ,4>{0.6830127f,1.1830127f,0.3169873f,-0.1830127f};
	constexpr static auto haahr = vector_td<typename realType<T>::Type,2>{1.0f,1.0f};
	constexpr static auto daubechies6= vector_td<typename realType<T>::Type,6>{0.47046721f,1.14111692f,0.650365f,-0.19093442f, -0.12083221f,0.0498175f};

	virtual boost::shared_ptr< linearOperator< cuNDArray<T> > > clone()
    								{
		return linearOperator<cuNDArray<T>>::clone(this);
    								}


private:

	std::vector<size_t> run_dimensions;
	unsigned int levels=0;
	bool use_random_ = false;


};
}
