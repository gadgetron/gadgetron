#pragma once

#include "hoNDArray_math.h"
#include "generalOperator.h"

#include "vector_td_operators.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron{

template<class T, unsigned int D> class hoTvOperator
		: public generalOperator< hoNDArray<T> >
{
protected:
	typedef typename realType<T>::Type REAL;

public:
	hoTvOperator() : generalOperator< hoNDArray<T> >(){
		limit_ = REAL(1e-8);
	}

	virtual ~hoTvOperator() {}

	void set_limit(REAL limit){
		limit_ = limit;
	}

	virtual void gradient( hoNDArray<T> *in_array, hoNDArray<T> *out_array, bool accumulate=false )
	{
		if (in_array->get_number_of_elements() != out_array->get_number_of_elements()){
			throw std::runtime_error("hoTvOperator: input/output array dimensions mismatch");
		}

		T* in = in_array->get_data_ptr();
		T* out = out_array->get_data_ptr();

		auto dims = vector_td<int64_t,D>(from_std_vector<size_t, D>(in_array->dimensions()));

		if (!accumulate)
			clear(out_array);

#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int64_t idx=0; idx < in_array->get_number_of_elements(); idx++){

			T xi = in[idx];
			T result = T(0);

			vector_td<int64_t,D> co = idx_to_co(idx, dims);

			REAL grad = gradient_(in,dims,co);

			if (grad > limit_) {
				result += REAL(D)*xi/grad;
				for (int i = 0; i < D; i++){
					co[i]+=1;
					result -= in[co_to_idx((co+dims)%dims,dims)]/grad;
					co[i]-=1;
				}
			}

			for (int i = 0; i < D; i++){
				co[i]-=1;
				grad = gradient_(in,dims,co);
				if (grad > limit_) {
					result +=(xi-in[co_to_idx((co+dims)%dims,dims)])/grad;
				}
				co[i]+=1;
			}
			out[idx] += this->weight_*result;
		}
	}


	virtual REAL magnitude( hoNDArray<T> *in_array )
	{

		T* in = in_array->get_data_ptr();

		 auto dims = vector_td<int64_t,D>(from_std_vector<size_t, D>(in_array->dimensions()));

		REAL result =0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:result)
#endif
		for (int64_t idx=0; idx < in_array->get_number_of_elements(); idx++){
			auto co = idx_to_co(idx, dims);
			REAL grad = gradient_(in,dims,co);
			result += this->weight_*grad;
		}

		return result;
	}

private:

	REAL inline gradient_(T* in, const vector_td<int64_t,D> dims, vector_td<int64_t,D> co)
	{
		REAL grad = REAL(0);
		T xi = in[co_to_idx((co+dims)%dims,dims)];
		for (int i = 0; i < D; i++){
			co[i]+=1;
			T dt = in[co_to_idx((co+dims)%dims,dims)];
			grad += norm(xi-dt);
			co[i]-=1;
		}
		return std::sqrt(grad);
	}

protected:
	REAL limit_;
};
}
