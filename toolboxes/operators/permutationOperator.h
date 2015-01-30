#pragma once

#include "linearOperator.h"
namespace Gadgetron{


template <class ARRAY_TYPE> class permutationOperator : public linearOperator<ARRAY_TYPE>{
public:

	virtual ~permutationOperator(){};

	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		if (accumulate)
			*out += *permute(in,&order);
		else
			permute(in,out,&order);
	}

	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		if (accumulate)
			*out += *permute(in,&transpose_order);
		else
			permute(in,out,&transpose_order);
	}

	virtual void mult_MH_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		if (accumulate)
			*out += *in;
		else
			*out = *in;
	}

	virtual void set_order(std::vector<size_t> order){
		this->order =order;
		transpose_order = std::vector<size_t>(order.size(),0);
		for (unsigned int i = 0; i < order.size(); i++)
			transpose_order[order[i]] = i;
	}


  virtual boost::shared_ptr< linearOperator<ARRAY_TYPE> > clone() {
    return linearOperator<ARRAY_TYPE>::clone(this);
  }
protected:
	 std::vector<size_t> order;
	 std::vector<size_t> transpose_order;
};
}
