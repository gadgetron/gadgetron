#pragma once

#include "cuTvOperator.h"
#include "hoCuNDArray.h"
#include "generalOperator.h"
#include "hoNDArray_fileio.h"
#include "vector_td_utilities.h"

namespace Gadgetron{

template<class T, unsigned int D> class EXPORTGPUOPERATORS hoCuTvOperator : public generalOperator<hoCuNDArray<T> >{
private:
	typedef typename realType<T>::Type REAL;
public:
	hoCuTvOperator() : generalOperator< hoCuNDArray<T> >(){
		limit_ = REAL(1e-8);
		cuTV.set_limit(limit_);
	}
	virtual ~hoCuTvOperator(){}
	void set_limit(REAL limit){
		limit_ = limit;
		cuTV.set_limit(limit);
	}

	virtual void gradient(hoCuNDArray<T>* in ,hoCuNDArray<T>* out, bool accumulate=false){
		if (in->get_number_of_elements() != out->get_number_of_elements()){
			BOOST_THROW_EXCEPTION(runtime_error("in_array and output array dimension mismatch"));
		}


		const vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in->get_dimensions()));
		int elements = in->get_number_of_elements();

		for (int i =0; i < (elements/prod(dims)); i++){
			std::vector<unsigned int> dimensions = to_std_vector(dims);
			hoNDArray<T> tmp_in;
			tmp_in.create(&dimensions,in->get_data_ptr()+i*prod(dims));
			hoNDArray<T> tmp_out;
			tmp_out.create(&dimensions,out->get_data_ptr()+i*prod(dims));
			cuNDArray<T> cuIn(&tmp_in);
			cuNDArray<T> cuOut(&tmp_out);

			cuTV.gradient(&cuIn,&cuOut,accumulate);
			boost::shared_ptr< hoNDArray<T> > tmp = cuOut.to_host();
			tmp_out = *tmp;
		}


	}

	virtual void set_weight(REAL weight){
		this->weight_ = weight;
		cuTV.set_weight(weight);
	}



protected:
	REAL limit_;
	cuTvOperator<T,D> cuTV;
};

template<class T> class hoCuTvOperator<T,4> : public generalOperator< hoCuNDArray<T> >{
private:
	typedef typename realType<T>::Type REAL;
public:
	hoCuTvOperator() : generalOperator< hoCuNDArray<T> >(){
		limit_ = REAL(1e-8);
	}
	virtual ~hoCuTvOperator(){}
	void set_limit(REAL limit){
		limit_ = limit;
	}

	virtual void gradient(hoCuNDArray<T>* in_array ,hoCuNDArray<T>* out_array, bool accumulate=false){
		if (in_array->get_number_of_elements() != out_array->get_number_of_elements()){
			BOOST_THROW_EXCEPTION(runtime_error("in_array and output array dimension mismatch"));
		}
		T* in = in_array->get_data_ptr();
		T* out = out_array->get_data_ptr();
		vector_td<unsigned int,4> dims = from_std_vector<unsigned int, 4>(*(in_array->get_dimensions()));
		if (!accumulate)clear(in_array);
		for (int idx = 0; idx < in_array->get_number_of_elements(); idx++){
				T xi = in[idx];
				T result=T(0);

				vector_td<unsigned int,4> co = idx_to_co<4>(idx, dims);

				REAL grad = gradient_(in,dims,co);


				if (grad > limit_) {
					result += REAL(4)*xi/grad;
					for (int i = 0; i < 4; i++){
							co[i]+=1;
							result -= in[co_to_idx<4>((co+dims)%dims,dims)]/grad;
							co[i]-=1;
						}
				}

				for (int i = 0; i < 4; i++){
					co[i]-=1;
					grad = gradient_(in,dims,co);
					if (grad > limit_) {
						result +=(xi-in[co_to_idx<4>((co+dims)%dims,dims)])/grad;
					}
					co[i]+=1;
				}
				out[idx] += this->weight_*result;

			}
	}

private:
	REAL inline gradient_(T* in, const vector_td<unsigned int,4> dims, vector_td<unsigned int,4> co){
		REAL grad = REAL(0);
		T xi = in[co_to_idx<4>((co+dims)%dims,dims)];
		for (int i = 0; i < 4; i++){
			co[i]+=1;
			T dt = in[co_to_idx<4>((co+dims)%dims,dims)];
			grad += norm(xi-dt);
			co[i]-=1;
		}
		return std::sqrt(grad);
	}


protected:
	REAL limit_;


};
}
