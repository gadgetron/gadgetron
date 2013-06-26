/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container)
*/

#ifndef HONDARRAY_H
#define HONDARRAY_H
#pragma once

#include "NDArray.h"
#include "complext.h"
#include "vector_td.h"

#include <string.h>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <new>          // std::runtime_error
namespace Gadgetron{

  template <class T> class hoNDArray : public NDArray<T>
  {
  public:

    hoNDArray() : NDArray<T>::NDArray() {}

    hoNDArray(std::vector<unsigned int> *dimensions) : NDArray<T>::NDArray() {
      this->create(dimensions);
    }

    hoNDArray(std::vector<unsigned int> *dimensions, 
	      T* data, bool delete_data_on_destruct = false) : NDArray<T>::NDArray() {
      this->create(dimensions,data,delete_data_on_destruct);
    }

    hoNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions) : NDArray<T>::NDArray() {
      this->create(dimensions.get());
    }

    hoNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions, 
	      T* data, bool delete_data_on_destruct = false) : NDArray<T>::NDArray() {
      this->create(dimensions.get(),data,delete_data_on_destruct);
    }
    
    virtual ~hoNDArray() {
      if (this->delete_data_on_destruct_) {
	deallocate_memory();
      }
    }

    // Copy constructor
    hoNDArray(const hoNDArray<T>& a) {
      this->data_ = 0;
      this->dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(new std::vector<unsigned int>(*a.dimensions_));
      allocate_memory();
      memcpy( this->data_, a.data_, this->elements_*sizeof(T) );
    }

    // Assignment operator
    hoNDArray& operator=(const hoNDArray& rhs) {
      // Are the dimensions the same? Then we can just memcpy
      if (this->dimensions_equal(&rhs)) {
	memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
      } else {
	deallocate_memory();
	this->data_ = 0;
	this->dimensions_ = rhs.dimensions_;
	allocate_memory();
	memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );

      }
      return *this;
    }
    
    virtual void fill(T value) {
#ifdef USE_OMP
#pragma omp parallel for
      for( int i=0; i<this->get_number_of_elements(); i++ )
	this->get_data_ptr()[i] = value;
#else
      std::fill(this->get_data_ptr(), this->get_data_ptr()+this->get_number_of_elements(), value);
#endif
    }
    
    T* begin() {
      return this->data_;
    }
    
    T* end() {
      return (this->data_+this->elements_);
    }
    
    T& at( unsigned int idx ){
      if( idx >= this->get_number_of_elements() ){
  	throw std::runtime_error("cuNDArray::at(): index out of range.");
      }
      return this->get_data_ptr()[idx];
    }
    
    T& operator[]( unsigned int idx ){
      if( idx >= this->get_number_of_elements() ){
  	throw std::runtime_error("cuNDArray::operator[]: index out of range.");
      }
      return this->get_data_ptr()[idx];
    }
    
  protected:
  
    virtual void allocate_memory()
    {
      deallocate_memory();

      this->elements_ = 1;
      for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
	this->elements_ *= (*this->dimensions_)[i];
      }

      _allocate_memory(this->elements_, &this->data_);
    
      if( this->data_ == 0x0 ){
	throw std::runtime_error("hoNDArray<>::allocate memory failed");
      }
    }

    virtual void deallocate_memory() {
      if( this->data_ ) {
	_deallocate_memory( this->data_ );
	this->data_ = 0x0;
      }
    }

    // Generic allocator / deallocator
    //

    template<class X> inline void _allocate_memory( unsigned int size, X** data ){
      *data = new (std::nothrow) X[size];
    }

    template<class X> inline void _deallocate_memory( X* data ){
      delete [] data;
    }
      
    // Overload these instances to avoid invoking the element class constructor/destructor
    //
    
    void inline _allocate_memory( unsigned int size, float_complext** data ){
      *data = (float_complext*) malloc( size*sizeof(float_complext) );
    }

    void inline _deallocate_memory( float_complext* data ){
      free( data );
    }

    void inline _allocate_memory( unsigned int size, double_complext** data ){
      *data = (double_complext*) malloc( size*sizeof(double_complext) );
    }

    void inline _deallocate_memory( double_complext* data ){
      free( data );
    }

    template<class TYPE, unsigned int D> void inline _allocate_memory( unsigned int size, vector_td<TYPE,D>** data ){
      *data = (vector_td<TYPE,D>*) malloc( size*sizeof(vector_td<TYPE,D>) );
    }

    template<class TYPE, unsigned int D>  void inline _deallocate_memory( vector_td<TYPE,D>* data ){
      free( data );
    }
  };
}

#endif //HONDARRAY_H
