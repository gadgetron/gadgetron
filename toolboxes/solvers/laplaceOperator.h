#pragma once

#include "matrixOperator.h"
#include "vector_td.h"

template <class REAL, unsigned int D, class ARRAY_TYPE> 
class laplaceOperator : public matrixOperator<REAL, ARRAY_TYPE>
{
  
public:
  
  laplaceOperator( ) : matrixOperator<REAL,ARRAY_TYPE>() { }
  virtual ~laplaceOperator() {}
    
  virtual int mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return compute_laplace( in, out, accumulate );    
  }
  
  virtual int mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return compute_laplace( in, out, accumulate );
  }

  virtual int mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {    
    ARRAY_TYPE tmp(*in); 

    if( !tmp.get_data_ptr() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : intermediate memory allocation failed" << std::endl;
      return -1;
    }
    int res = mult_M(in,&tmp);

    if (res < 0) {
      return res;
    }
    
    return mult_MH(&tmp,out,accumulate);
  }
    
protected:
  virtual int compute_laplace( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;    
};
