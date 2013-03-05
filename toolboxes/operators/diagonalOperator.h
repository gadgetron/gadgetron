#pragma once

#include "linearOperator.h"
#include "vector_td_utilities.h"

#include <boost/smart_ptr.hpp>
#include <vector>

template <class REAL, class ARRAY_TYPE> class diagonalOperator 
  : public linearOperator<REAL, ARRAY_TYPE>
{
  
public:
  
  diagonalOperator() : linearOperator<REAL, ARRAY_TYPE>() {}
  virtual ~diagonalOperator() {}
  
  // Set/get diagonal
  virtual void set_diagonal( boost::shared_ptr<ARRAY_TYPE> diagonal ) { diagonal_ = diagonal; }
  virtual ARRAY_TYPE* get_diagonal() { return diagonal_.get(); }
      
  // Apply diagonal operator (twice)
  virtual int mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {    
    bool ret2 = true;
    
    if( !accumulate ) 
      ret2 = operator_clear( out );
    
    if( ret2 ){
      ret2 = operator_axpy( diagonal_.get(), in, out );
      ret2 &= operator_axpy( diagonal_.get(), in, out );
    }else
      ret2 = false;
        
    if( ret2 )
      return 0;
    else{
      std::cout << std::endl << "Error: diagonalOperator::mult_MH_M failed" << std::endl;
      return -1;
    }
  }
  
  virtual int mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    bool ret2 = true;

    if( !accumulate )
      ret2 = operator_clear( out );

    if( ret2 ){
      ret2 = operator_axpy( diagonal_.get(), in, out );

    }else
      ret2 = false;

    if( ret2 )
      return 0;
    else{
      std::cout << std::endl << "Error: diagonalOperator::mult_M failed" << std::endl;
      return -1;
    }
  }
  
  virtual int mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    bool ret2 = true;

    if( !accumulate )
      ret2 = operator_clear( out );

    if( ret2 ){
      ret2 = operator_axpy( diagonal_.get(), in, out );

    }else
      ret2 = false;

    if( ret2 )
      return 0;
    else{
      std::cout << std::endl << "Error: diagonalOperator::mult_MH failed" << std::endl;
      return -1;
    }
  }
  
  virtual bool operator_clear( ARRAY_TYPE* ) = 0;
  virtual bool operator_axpy( ARRAY_TYPE*, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  
protected:
  boost::shared_ptr<ARRAY_TYPE> diagonal_;
};
