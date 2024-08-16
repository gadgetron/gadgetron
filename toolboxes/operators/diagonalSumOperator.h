/** \file diagonalSumOperator.h
    \brief Operator to compute the sum over a set of diagonal matrices times a set of corresponding vectors.

    The domain of this operator is a set of images, the codomain a single image. 
    The sum is computed over the last dimension of the provided diagonal array.
*/

#pragma once

#include "diagonalOperator.h"

namespace Gadgetron {

  template <class ARRAY_TYPE> class diagonalSumOperator : public diagonalOperator<ARRAY_TYPE>
  {
  public:
  
    diagonalSumOperator() : diagonalOperator<ARRAY_TYPE>() {}
    virtual ~diagonalSumOperator() {}
  
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( !this->diagonal_ ){
        throw std::runtime_error("diagonalSumOperator::mult_M failed: diagonal not set");
      }       

      const unsigned int num_phases = this->diagonal_->get_size(this->diagonal_->get_number_of_dimensions()-1);
      const unsigned int elements_per_phase = this->diagonal_->get_number_of_elements()/num_phases;
      
      if( in->get_number_of_elements() != this->diagonal_->get_number_of_elements() ){
        throw std::runtime_error("diagonalSumOperator::mult_M failed: array size mismatch between input image and diagonal");
      }

      if( out->get_number_of_elements() != elements_per_phase ){
        throw std::runtime_error("diagonalSumOperator::mult_M failed: the output image domain should only be a single image");
      }

      if( !accumulate ) 
        clear(out);

      std::vector<size_t> dims = out->get_dimensions();
     
      // Iterate over the last dimension of the provided diagonal image
      //

      for( unsigned int i=0; i<num_phases; i++ ){

        ARRAY_TYPE tmp_in(dims, in->get_data_ptr()+i*elements_per_phase );
        ARRAY_TYPE tmp_diag(dims, this->diagonal_->get_data_ptr()+i*elements_per_phase );

        if(i==0 && !accumulate){
          *out = tmp_in;
          *out *= tmp_diag;
        }
        else{
          ARRAY_TYPE tmp(tmp_in);
          tmp *= tmp_diag;
          *out += tmp;
        }
      }
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( !this->diagonal_conj_ ){
        throw std::runtime_error("diagonalSumOperator::mult_MH failed: diagonal not set");
      }       

      const unsigned int num_phases = this->diagonal_conj_->get_size(this->diagonal_conj_->get_number_of_dimensions()-1);
      const unsigned int elements_per_phase = this->diagonal_conj_->get_number_of_elements()/num_phases;
      
      if( in->get_number_of_elements() != elements_per_phase ){
        throw std::runtime_error("diagonalSumOperator::mult_MH failed: the input image domain should only be a single image");
      }

      if( out->get_number_of_elements() != this->diagonal_conj_->get_number_of_elements() ){
        throw std::runtime_error("diagonalSumOperator::mult_MH failed: array size mismatch between output image and diagonal");
      }

      if( !accumulate ){
        *out = *this->diagonal_conj_;
        *out *= *in; // multiplies all phases with the input
      }
      else{
        ARRAY_TYPE tmp(*this->diagonal_conj_);
        tmp *= *in; // multiplies all phases with the input
        *out += tmp;
      }
    }

  };
}
