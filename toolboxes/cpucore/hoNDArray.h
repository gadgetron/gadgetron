#ifndef HONDARRAY_H
#define HONDARRAY_H
#pragma once

#include "NDArray.h"

#include "complext.h"
#include "vector_td.h"

#include <string.h>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

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

template <class T> class hoNDArray : public NDArray<T>
{
public:

  hoNDArray() : NDArray<T>::NDArray() {}

  hoNDArray(std::vector<unsigned int> *dimensions) : NDArray<T>::NDArray() {
	  this->create(dimensions);
  }
  hoNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray() {
  	  this->create(dimensions,data,delete_data_on_destruct);
    }

  hoNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions) : NDArray<T>::NDArray() {
	  this->create(dimensions.get());
  }
  hoNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray() {
  	  this->create(dimensions.get(),data,delete_data_on_destruct);
    }
  virtual ~hoNDArray() {
    if (this->delete_data_on_destruct_) {
      deallocate_memory();
    }
  }



  // Copy constructor
  hoNDArray(const hoNDArray<T>& a)
  {
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

  static boost::shared_ptr< hoNDArray<T> > allocate(std::vector<unsigned int> *dimensions)
  {
    boost::shared_ptr< hoNDArray<T> > ret( new hoNDArray<T> );

    if( ret->create(dimensions) == 0x0 ) {
      std::cerr << "hoNDArray<T>::allocate failed to create array" << std::endl;
      return boost::shared_ptr< hoNDArray<T> >();
    }

    return ret;
  }

  virtual void clear(T value)
  {
    std::fill(this->get_data_ptr(), this->get_data_ptr()+this->get_number_of_elements(), value);
  }

  virtual void permute(std::vector<unsigned int> *dim_order, NDArray<T>* out = 0, int shift_mode = 0)
  {
    hoNDArray<T>* out_int = 0;

    // Check ordering array
    if (dim_order->size() > this->dimensions_->size()) {
      std::runtime_error("hoNDArray::permute - Invalid length of dimension ordering array");
    }

    std::vector<unsigned int> dim_count(this->dimensions_->size(),0);
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if ((*dim_order)[i] >= this->dimensions_->size()) {
	std::runtime_error("hoNDArray::permute - Invalid dimension order array");

      }
      dim_count[(*dim_order)[i]]++;
    }

    // Create an internal array to store the dimensions
    std::vector<unsigned int> dim_order_int;

    // Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if (dim_count[(*dim_order)[i]] != 1) {
	std::runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)");

      }
      dim_order_int.push_back((*dim_order)[i]);
    }

    // Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < this->dimensions_->size()) {
      for (unsigned int i = 0; i < dim_count.size(); i++) {
	if (dim_count[i] == 0) {
	  dim_order_int.push_back(i);
	}
      }
    }

    if (out) {
      out_int = dynamic_cast< hoNDArray<T>* >(out);
      if (!out_int) {
    	  std::runtime_error("hoNDArray::permute: failed to dynamic cast out array pointer");

      }
      for (unsigned int i = 0; i < dim_order_int.size(); i++) {
	if ((*this->dimensions_)[dim_order_int[i]] != out_int->get_size(i)) {
	  std::runtime_error("hoNDArray::permute: Dimensions of output array do not match the input array");
	}
      }
    }

    T* o = 0;
    if (out_int) {
      o = out_int->get_data_ptr();
    } else {
      o = new T[this->elements_];
    }

    ArrayIterator it(this->dimensions_.get(),&dim_order_int);
    for (unsigned long int i = 0; i < this->elements_; i++) {
      o[i] = this->data_[it.get_current_idx()];
      it.advance();
    }

    if (!out_int) {
      boost::shared_ptr<std::vector<unsigned int> > tmp_dims(new std::vector<unsigned int>);
      for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
	tmp_dims->push_back((*this->dimensions_)[dim_order_int[i]]);
      }
      this->dimensions_ = tmp_dims;

      delete [] this->data_;
      this->data_ = o;
    }


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
      std::runtime_error("hoNDArray<>::allocate memory failed");
    }
    
  }

  virtual void deallocate_memory() {
    if( this->data_ ) {
      _deallocate_memory( this->data_ );
      this->data_ = 0x0;
    }

  }

  // Generic allocator / deallocator
  template<class X> inline void _allocate_memory( unsigned int size, X** data ){
    *data = new X[size];
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

#endif //HONDARRAY_H
