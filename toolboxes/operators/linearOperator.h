/** \file linearOperator.h
    \brief Base class for all linear operators.
*/

#pragma once

#include "generalOperator.h"

namespace Gadgetron{

  /** \class linearOperator
      \brief Base class for all linear Operators
  */
  template <class ARRAY_TYPE> class linearOperator : public generalOperator<ARRAY_TYPE>
  {
  public:
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
  linearOperator() : generalOperator<ARRAY_TYPE>() {}

  linearOperator(std::vector<size_t> *dims) : generalOperator<ARRAY_TYPE>(dims) {
      set_codomain_dimensions(dims);
    }

  linearOperator(std::vector<size_t> *dims, std::vector<size_t> *codims)
    : generalOperator<ARRAY_TYPE>(dims) {
      set_codomain_dimensions(codims);
    }

    virtual ~linearOperator() {}

    /**
     * The gradient of a linear operator corresponds to mult_MH_M, times the weight of the operator.
     * @param[in] in Input array.
     * @param[in,out] out Output Array.
     * @param accumulate If true, adds result to out. If false, overwrites out.
     */
    virtual void gradient(ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false)
    {
      if( in == 0x0 || out == 0x0 )
	throw std::runtime_error("linearOperator::gradient(): Invalid input and/or output array");

      ARRAY_TYPE* tmp = out;
      if (accumulate) {
	tmp = new ARRAY_TYPE(out->get_dimensions());
      }
      mult_MH_M(in,tmp,false);
      *tmp *= this->weight_;
      if (accumulate){
	*out += *tmp;
	delete tmp;
      }
    }


    virtual REAL magnitude(ARRAY_TYPE* in){
      ARRAY_TYPE tmp(&this->codomain_dims_);
      this->mult_M(in,&tmp);
      return std::sqrt(this->get_weight())*real(dot(&tmp,&tmp));
    }
    virtual void set_codomain_dimensions( const std::vector<size_t> *dims )
    {
      if( dims == 0x0 )
	throw std::runtime_error("linearOperator::set_codomain_dimensions: illegal dimensions array provided");
      codomain_dims_ = *dims;
    }

    virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions()
      {
	std::vector<size_t> *dims = new std::vector<size_t>();
	*dims = codomain_dims_;
	return boost::shared_ptr< std::vector<size_t> >(dims);
      }

    virtual void mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
    virtual void mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;

    virtual void mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
    {
      if( codomain_dims_.size() == 0 ){
	throw std::runtime_error("Error: linearOperator::mult_MH_M : codomain dimensions not set");
      }

      ARRAY_TYPE tmp;
      tmp.create(&codomain_dims_);
      mult_M( in, &tmp, false );
      mult_MH( &tmp, out, accumulate );
    }

  protected:
    std::vector<size_t> codomain_dims_;
  };
}
