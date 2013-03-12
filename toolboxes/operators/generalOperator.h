#pragma once
#include "GadgetronException.h"
namespace Gadgetron{
/**
 * A base class for any sort of optimization parameter on which we can calculate a gradient.
 */
template <class ARRAY> class generalOperator{
private:
		typedef typename ARRAY::element_type ELEMENT_TYPE;
	  typedef typename realType<ELEMENT_TYPE>::type REAL;
public:
	generalOperator():weight_(REAL(1)){};
	generalOperator(std::vector<unsigned int> *dims):weight_(REAL(1)){
		set_domain_dimensions(dims);
	}
	virtual void gradient(ARRAY* in, ARRAY* out,bool accumulate=false)=0;
	virtual ~generalOperator(){};
	virtual void set_domain_dimensions( std::vector<unsigned int> *dims )
	{
		if( dims == 0x0 ) BOOST_THROW_EXCEPTION(runtime_error("Null pointer provided"));
		domain_dims_ = *dims;

	}
 virtual void set_weight( REAL weight ){ weight_ = weight; }
 virtual REAL get_weight(){ return weight_; }
 virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions()
 {
   std::vector<unsigned int> *dims = new std::vector<unsigned int>();
   *dims = domain_dims_;
   return boost::shared_ptr< std::vector<unsigned int> >(dims);
 }
protected:
 REAL weight_;
	std::vector<unsigned int> domain_dims_;
};

}
