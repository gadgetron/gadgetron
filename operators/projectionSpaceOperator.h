#pragma once
#include "linearOperator.h"
#include <boost/shared_ptr.hpp>

namespace Gadgetron{
template<class ARRAY_TYPE> class projectionSpaceOperator : public linearOperator<ARRAY_TYPE>{
public:
	projectionSpaceOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE > > _op) : linearOperator<ARRAY_TYPE>(), op(_op){};
	projectionSpaceOperator() : linearOperator<ARRAY_TYPE>() {};
	void set_projections(boost::shared_ptr<ARRAY_TYPE > _projections){projections = _projections;}
	virtual ~projectionSpaceOperator(){};
	virtual void mult_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		op->mult_M(in,out,accumulate);
	}
	virtual void mult_MH(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		op->mult_M(in,out,accumulate);
	}

	virtual void mult_MH_M(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate=false){
		op->mult_M(in,out,accumulate);
	}

	 /**
		 * The gradient of a linear operator corresponds to mult_MH_M.
		 * @param[in] in Input array.
		 * @param[in,out] out Output Array.
		 * @param accumulate If true, adds result to out. If false, overwrites out.
		 */
	virtual void gradient(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false)
	{
		if( in == 0x0 || out == 0x0 )
			BOOST_THROW_EXCEPTION(runtime_error("linearOperator::gradient(): Invalid input and/or output array"));

		ARRAY_TYPE* tmp = out;
		if (accumulate) {
			tmp = new ARRAY_TYPE(out->get_dimensions());
		}

		ARRAY_TYPE pSpace(projections->get_dimensions());

		this->mult_M(in,&pSpace,false);
		pSpace -= *projections;
		this->mult_MH(&pSpace,tmp,false);

		*tmp *= this->weight_;
		if (accumulate){
			*out += *tmp;
			delete tmp;
		}
	}

		virtual boost::shared_ptr< linearOperator<ARRAY_TYPE > > clone(){
			return linearOperator<ARRAY_TYPE>::clone(this);
		}
protected:

	boost::shared_ptr<linearOperator<ARRAY_TYPE > > op;

	boost::shared_ptr<ARRAY_TYPE > projections;
};
}
