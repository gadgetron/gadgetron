#pragma once

#include "linearOperator.h"
namespace Gadgetron{

template<class ARRAY_TYPE> class weightingOperator: public linearOperator<ARRAY_TYPE>{
public:
	weightingOperator(): linearOperator<ARRAY_TYPE>(){};
	weightingOperator(boost::shared_ptr<ARRAY_TYPE> _weights): linearOperator<ARRAY_TYPE>(), weights(_weights){};

	virtual void set_weights(boost::shared_ptr<ARRAY_TYPE> _weights){
		weights = _weights;
	}

	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		if (accumulate){
			ARRAY_TYPE tmp = *in;
			tmp *= *weights;
			*out += tmp;
		} else{
			*out = *in;
			*out *= *weights;
		}
	}

	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		mult_M(in,out,accumulate);
	}

	virtual void mult_MH_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
			if (accumulate){
				ARRAY_TYPE tmp = *in;
				tmp *= *weights;
				tmp *= *weights;
				*out += tmp;
			} else{
				*out = *in;
				*out *= *weights;
				*out *= *weights;
			}
		}
	virtual boost::shared_ptr< linearOperator<ARRAY_TYPE > > clone(){
				return linearOperator<ARRAY_TYPE>::clone(this);
			}
protected:
	boost::shared_ptr<ARRAY_TYPE> weights;
};

}
