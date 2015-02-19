/** \file diagonalOperator.h
    \brief Base class for the diagonal matrix operators.
*/

#pragma once

#include "linearOperator.h"

namespace Gadgetron {

  template <class ARRAY_TYPE> class diagonalOperator : public linearOperator<ARRAY_TYPE>
  {
  public:
  
    diagonalOperator() : linearOperator<ARRAY_TYPE>() {}
    virtual ~diagonalOperator() {}
  
    // Set/get diagonal
    //
    
    virtual void set_diagonal( boost::shared_ptr<ARRAY_TYPE> diagonal ) { 
      diagonal_ = diagonal;
      diagonal_conj_ = conj(diagonal.get());
    }

    virtual boost::shared_ptr<ARRAY_TYPE> get_diagonal() { return diagonal_; }
  
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ) {
        ARRAY_TYPE tmp(*in);
        tmp *= *diagonal_;
        *out += tmp;
      }
      else{
        *out = *in;
        *out *= *diagonal_;
      }
    }
  
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ) {
        ARRAY_TYPE tmp(*in);
        tmp *= *diagonal_conj_;
        *out += tmp;
      }
      else{
        *out = *in;
        *out *= *diagonal_conj_;
      }
    }
    
    // Apply diagonal operator (twice)
    virtual void mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate )
    {    
      if( accumulate ) {
        ARRAY_TYPE tmp(*in);
        tmp *= *diagonal_;
        tmp *= *diagonal_conj_;
        *out += tmp;
      }
      else{
        *out = *in;
        *out *= *diagonal_;
        *out *= *diagonal_conj_;
      }
    }
  
  protected:
    boost::shared_ptr<ARRAY_TYPE> diagonal_;
    boost::shared_ptr<ARRAY_TYPE> diagonal_conj_;
  };
}
