#pragma once
#include "linearOperator.h"
#include <numeric>
#include <functional>
namespace Gadgetron {

  template<class ARRAY_TYPE,class partialDerivOp, unsigned int D> class opticalFlowOperator : public linearOperator<ARRAY_TYPE >{
  public:
    typedef typename ARRAY_TYPE::element_type T;

    opticalFlowOperator(){};

    opticalFlowOperator(ARRAY_TYPE* moving,ARRAY_TYPE* stat){
      set_images(moving,stat);
    }

    virtual ~opticalFlowOperator(){};

    virtual void mult_M(ARRAY_TYPE* in,ARRAY_TYPE* out,bool accumulate){

      if (!accumulate) clear(out);
      std::vector<size_t> dims = in->get_dimensions();
      if (dims.back() != D) throw std::runtime_error("Input array for optical flow has the wrong last dimensions");
      dims.pop_back();

      size_t elements = std::accumulate(dims.begin(),dims.end(),1u,std::multiplies<size_t>());

      for (int i = 0; i < D; i++){
	ARRAY_TYPE tmp(dims,in->get_data_ptr()+elements*i);
	ARRAY_TYPE tmp2(tmp);
	tmp2 *= *Ix[i];
	*out += tmp2;
      }
    }

    virtual void mult_MH(ARRAY_TYPE* in,ARRAY_TYPE* out,bool accumulate){

      if (!accumulate) clear(out);
      std::vector<size_t> dims = out->get_dimensions();
      if (dims.back() != D) throw std::runtime_error("Output array for optical flow has the wrong last dimensions");
      dims.pop_back();
      size_t elements = std::accumulate(dims.begin(),dims.end(),1u,std::multiplies<size_t>());

      for (int i = 0; i < D; i++){
	ARRAY_TYPE out_view(dims,out->get_data_ptr()+elements*i);
	ARRAY_TYPE tmp2(*in);
	tmp2 *= *Ix[i];
	out_view += tmp2;
      }
    }

    void set_images(ARRAY_TYPE* moving,ARRAY_TYPE* stat){
      Ix = std::vector< boost::shared_ptr<ARRAY_TYPE> >();

      for (int i=0; i < D; i++){
	partialDerivOp op(i);
	boost::shared_ptr<ARRAY_TYPE> I(new ARRAY_TYPE(moving->get_dimensions()));
	op.mult_M(moving,I.get());
	op.mult_M(stat,I.get(),true);
	*I /= T(2);
	Ix.push_back(I);
      }
    }

  protected:
    std::vector< boost::shared_ptr<ARRAY_TYPE> > Ix; //Gradient along different directions
  };
}
