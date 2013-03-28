/** \file NDArray.h
    \brief Abstract base class for all Gadgetron host and device arrays
*/

#ifndef NDARRAY_H
#define NDARRAY_H
#pragma once

#include "GadgetronException.h"

#include <new>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  template <class T> class NDArray
  {
  public:
    typedef T element_type;

    void* operator new (size_t bytes) {
      return ::new char[bytes];
    }
  
    void operator delete (void *ptr) {
      delete [] static_cast <char *> (ptr);
    } 
  
    void * operator new(size_t s, void * p) {
      return p;
    }
  
    NDArray () : data_(0), elements_(0), delete_data_on_destruct_(true) {
      dimensions_ = boost::shared_ptr< std::vector<unsigned int> >( new std::vector<unsigned int> );
    }
  
    virtual ~NDArray() {}
  
    virtual void create(std::vector<unsigned int> *dimensions) 
    {
      std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
      *tmp = *dimensions;
      dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(tmp);
      allocate_memory();
    }
    
    virtual void create(std::vector<unsigned int> *dimensions, 
			T* data, bool delete_data_on_destruct = false) 
    {
      if (!data) {
	BOOST_THROW_EXCEPTION(runtime_error("NDArray<T>::create: 0x0 pointer provided"));
      }
      
      std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
      *tmp = *dimensions;
      dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(tmp);
      this->data_ = data;
      this->delete_data_on_destruct_ = delete_data_on_destruct;
      this->elements_ = 1;
      for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
	this->elements_ *= (*this->dimensions_)[i];
      }
    }
    
    virtual void create(boost::shared_ptr< std::vector<unsigned int> > dimensions) {
      this->create(dimensions.get());
    }
  
    virtual void create(boost::shared_ptr<std::vector<unsigned int>  > dimensions, 
			T* data, bool delete_data_on_destruct = false) {
      this->create(dimensions.get(), data, delete_data_on_destruct);
    }
      
    inline void squeeze() {
      boost::shared_ptr< std::vector<unsigned int> > new_dimensions( new std::vector<unsigned int> ); 
      for (unsigned int i = 0; i < dimensions_->size(); i++) {
	if ((*dimensions_)[i] != 1) {
	  new_dimensions->push_back((*dimensions_)[i]);
	}
      }
      dimensions_ = new_dimensions;
    }
  
    inline void reshape(std::vector<unsigned int> *dims) {
      unsigned long int new_elements = 1;
      for (unsigned int i = 0; i < dims->size(); i++) {
	new_elements *= (*dims)[i];
      }
    
      if (new_elements != elements_) {
	std::runtime_error("NDArray<T>::reshape : Number of elements cannot change during reshape");
      }
    
      // Copy the input dimensions array
      std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
      *tmp = *dims;
      dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(tmp);
    }
  
    inline void reshape( boost::shared_ptr< std::vector<unsigned int> > dims ) {
      reshape(dims.get());
    }

    inline bool dimensions_equal(std::vector<unsigned int> *d) {
      return ((this->dimensions_->size() == d->size()) &&
	      std::equal(this->dimensions_->begin(), this->dimensions_->end(), d->begin()));
    }
    
    template<class S> bool dimensions_equal(const NDArray<S> *a) {
      boost::shared_ptr<std::vector<unsigned int > > adims = a->get_dimensions();
      return ((this->dimensions_->size() == adims->size()) &&
	      std::equal(this->dimensions_->begin(), this->dimensions_->end(), adims->begin()));
    }
  
    inline unsigned int get_number_of_dimensions() const {
      return dimensions_->size();
    }
  
    inline unsigned int get_size(unsigned int dimension) const {
      if (dimension >= dimensions_->size()) {
	return 1;
      } else {
	return (*dimensions_)[dimension];
      }
    }
  
    inline boost::shared_ptr< std::vector<unsigned int> > get_dimensions() const {
      // Make copy to ensure that the receiver cannot alter the array dimensions
      std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
      *tmp=*dimensions_;
      return boost::shared_ptr< std::vector<unsigned int> >(tmp); 
    }
  
    inline T* get_data_ptr() const { 
      return data_; 
    }
  
    inline unsigned long int get_number_of_elements() const {
      return elements_;
    }
  
    inline bool delete_data_on_destruct() const {
      return delete_data_on_destruct_;
    }
  
    inline void delete_data_on_destruct(bool d) {
      delete_data_on_destruct_ = d;
    }

  protected:

    virtual void allocate_memory() = 0;
    virtual void deallocate_memory() = 0;

  protected:

    boost::shared_ptr< std::vector<unsigned int> > dimensions_;
    T* data_;
    unsigned long int elements_;
    bool delete_data_on_destruct_;  
  };
}

#endif //NDARRAY_H
