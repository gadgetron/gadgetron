/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container)
*/

#ifndef HONDARRAY_H
#define HONDARRAY_H
#pragma once

#include "NDArray.h"
#include "complext.h"
#include "vector_td.h"
#include "GadgetronException.h"

#include <string.h>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

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
      this->dimensions_ = a.dimensions_;
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

    class ArrayIterator
    {
    public:
      
      ArrayIterator(std::vector<unsigned int> *dimensions, std::vector<unsigned int> *order)
      {
	dimensions_  = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
	order_       = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
	current_     = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
	block_sizes_ = boost::shared_ptr< std::vector<unsigned long int> > (new std::vector<unsigned long int>);
	
	block_sizes_->push_back(1);
	for (unsigned int i = 0; i < order->size(); i++) {
	  dimensions_->push_back((*dimensions)[i]);
	  order_->push_back((*order)[i]);
	  current_->push_back(0);
	  if (i > 0) {
	    block_sizes_->push_back((*block_sizes_)[i-1]*(*dimensions_)[i-1]);
	  }
	}
	current_idx_ = 0;
      }
      
      inline unsigned long int advance()
      {
	unsigned int order_index = 0;
	(*current_)[(*order_)[order_index]]++;
	while ((*current_)[(*order_)[order_index]] >= (*dimensions_)[(*order_)[order_index]]) {
	  (*current_)[(*order_)[order_index]] = 0;
	  order_index = (order_index+1)%dimensions_->size();
	  (*current_)[(*order_)[order_index]]++;
	}
	
	current_idx_ = 0;
	for (unsigned int i = 0; i < dimensions_->size(); i++) {
	  current_idx_ += (*current_)[i]*(*block_sizes_)[i];
	}
	
	return current_idx_;
      }
      
      inline unsigned long int get_current_idx()
      {
	return current_idx_;
      }
      
      boost::shared_ptr< std::vector<unsigned int> > get_current_sub()
      {
	return current_;
      }
      
    protected:
      boost::shared_ptr< std::vector<unsigned int> > dimensions_;
      boost::shared_ptr< std::vector<unsigned int> > order_;
      boost::shared_ptr< std::vector<unsigned int> > current_;
      boost::shared_ptr< std::vector<unsigned long int> > block_sizes_;
      unsigned long int current_idx_;
    };

    T* begin() {
      return this->data_;
    }
    
    T* end() {
      return (this->data_+this->elements_);
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
	BOOST_THROW_EXCEPTION( bad_alloc("hoNDArray<>::allocate memory failed"));
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
