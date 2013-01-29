#ifndef NDARRAY_H
#define NDARRAY_H
#pragma once

#include <new>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace GADGETRON
{

template <class T> class NDArray
{
 public:

  void* operator new (size_t bytes)
  {
    return ::new char[bytes];
  }

  void operator delete (void *ptr)
  {
    delete [] static_cast <char *> (ptr);
  } 

  void * operator new(size_t s, void * p)
  {
    return p;
  }
  
  NDArray () : data_(0), elements_(0), delete_data_on_destruct_(true){
    dimensions_ = boost::shared_ptr< std::vector<unsigned int> >( new std::vector<unsigned int> );
  }
  
  virtual ~NDArray() {}
  
  virtual T* create(std::vector<unsigned int> *dimensions) 
  {
    std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
    *tmp = *dimensions;
    dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(tmp);
    allocate_memory();
    return this->get_data_ptr();
  }

  virtual T* create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false) 
  {
    if (!data) {
      std::cerr << "NDArray<T>::create: 0x0 pointer provided" << std::endl;
      return 0x0;
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
    
    return this->get_data_ptr();
  }

  virtual int permute(std::vector<unsigned int> *dim_order, NDArray<T> *out = 0, int shift_mode = 0) = 0;
  
  inline int shift_dim(int shift, NDArray<T> *out = 0) {
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < dimensions_->size(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%dimensions_->size()));
    }
    return permute(&order,out);
  }

  inline void squeeze(){
    boost::shared_ptr< std::vector<unsigned int> > new_dimensions( new std::vector<unsigned int> ); 
    for (unsigned int i = 0; i < dimensions_->size(); i++) {
      if ((*dimensions_)[i] != 1) {
	new_dimensions->push_back((*dimensions_)[i]);
      }
    }
    dimensions_ = new_dimensions;
  }

  inline int reshape(std::vector<unsigned int> *dims) {
    unsigned long int new_elements = 1;
    for (unsigned int i = 0; i < dims->size(); i++) {
      new_elements *= (*dims)[i];
    }
    
    if (new_elements != elements_) {
      std::cerr << "NDArray<T>::reshape : Number of elements cannot change during reshape" << std::endl;
      return -1;
    }

    // Copy the input dimensions array
    std::vector<unsigned int> *tmp = new std::vector<unsigned int>;
    *tmp = *dims;
    dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(tmp);
    return 0;
  }

  inline bool dimensions_equal(std::vector<unsigned int> *d) {
    return ((this->dimensions_->size() == d->size()) &&
	    std::equal(this->dimensions_->begin(), this->dimensions_->end(), d->begin()));
  }
  
  inline bool dimensions_equal(const NDArray<T> *a) {
    return ((this->dimensions_->size() == a->dimensions_->size()) &&
	    std::equal(this->dimensions_->begin(), this->dimensions_->end(), a->dimensions_->begin()));
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

  virtual int allocate_memory() = 0;
  virtual int deallocate_memory() = 0;

  boost::shared_ptr< std::vector<unsigned int> > dimensions_;
  T* data_;
  unsigned long int elements_;
  bool delete_data_on_destruct_;  
};

}

#endif //NDARRAY_H
