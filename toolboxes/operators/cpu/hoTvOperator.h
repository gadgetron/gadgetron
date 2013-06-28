#pragma once

#include "generalOperator.h"
#include "hoNDArray_operators.h"
#include "vector_td_operators.h"

#include <omp.h>

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

      vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in_array->get_dimensions()));

      if (!accumulate)
	clear(in_array);

#ifdef USE_OMP
#pragma omp parallel for
#endif
      for (int idx=0; idx < in_array->get_number_of_elements(); idx++){

	T xi = in[idx];
	T result = T(0);

	vector_td<unsigned int,D> co = idx_to_co<D>(idx, dims);

	REAL grad = gradient_(in,dims,co);

	if (grad > limit_) {
	  result += REAL(D)*xi/grad;
	  for (int i = 0; i < D; i++){
	    co[i]+=1;
	    result -= in[co_to_idx<D>((co+dims)%dims,dims)]/grad;
	    co[i]-=1;
	  }
	}

	for (int i = 0; i < D; i++){
	  co[i]-=1;
	  grad = gradient_(in,dims,co);
	  if (grad > limit_) {
	    result +=(xi-in[co_to_idx<D>((co+dims)%dims,dims)])/grad;
	  }
	  co[i]+=1;
	}
	out[idx] += this->weight_*result;
      }
    }
    
  private:

    REAL inline gradient_(T* in, const vector_td<unsigned int,D> dims, vector_td<unsigned int,D> co)
    {
      REAL grad = REAL(0);
      T xi = in[co_to_idx<D>((co+dims)%dims,dims)];
      for (int i = 0; i < D; i++){
	co[i]+=1;
	T dt = in[co_to_idx<D>((co+dims)%dims,dims)];
	grad += norm(xi-dt);
	co[i]-=1;
      }
      return std::sqrt(grad);
    }
    
  protected:
    REAL limit_;
  };
}
