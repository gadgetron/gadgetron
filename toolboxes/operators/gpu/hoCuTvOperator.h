#pragma once

#include "generalOperator.h"
#include "hoCuNDArray.h"
#include "cuTvOperator.h"
#include "vector_td_utilities.h"

namespace Gadgetron{

  template<class T, unsigned int D> class EXPORTGPUOPERATORS hoCuTvOperator : 
    public generalOperator< hoCuNDArray<T> >
  {

  protected:
    typedef typename realType<T>::Type REAL;

  public:
    
    hoCuTvOperator() : generalOperator< hoCuNDArray<T> >(){
      limit_ = REAL(1e-8);
      cuTV.set_limit(limit_);
    }

    virtual ~hoCuTvOperator(){}

    void set_limit( REAL limit ){
      limit_ = limit;
      cuTV.set_limit(limit);
    }
    
    virtual void gradient( hoCuNDArray<T> *in, hoCuNDArray<T> *out, bool accumulate=false )
    {
      if (in->get_number_of_elements() != out->get_number_of_elements()){
	BOOST_THROW_EXCEPTION(runtime_error("hoCuTvOperator: input/output array dimensions mismatch"));
      }
      
      const vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in->get_dimensions()));
      int elements = in->get_number_of_elements();
      
      for (int i=0; i < (elements/prod(dims)); i++){

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
}
