#ifndef NDARRAY_H
#define NDARRAY_H

#include <new>
#include <vector>
#include <iostream>

template <class T> class NDArray
{
 public:
  NDArray () 
    : data_(0)
    , elements_(0)
  {
    
  }

  ~NDArray() {
    deallocate_memory();
  }

  T* create(std::vector<int>& dimensions) {
    dimensions_ = dimensions; 
    allocate_memory();
    return get_data_ptr();
  }

  //Copy constructor
  NDArray(const NDArray<T>& a) {
    data_ = 0;
    dimensions_ = a.dimensions_;
    if (allocate_memory() == 0) {
      memcpy( data_, a.data_, elements_*sizeof(T) );
    }
  }
  
  T& operator[] (long int i);
  NDArray<T>& operator=(const NDArray<T> &a);
  
  int get_number_of_dimensions() {
    return dimensions_.size();
  }

  int get_size(unsigned int dimension) {
    if (dimension >= dimensions_.size() || dimension < 0) {
      return 1;
    } else {
      return dimensions_[dimension];
    }
  }

  std::vector<int> get_dimensions() {
    return dimensions_;
  }

  unsigned long int get_number_of_elements() {
    return elements_;
  }
  
  int clear() {
    if (data_ == 0) {
      return -1;
    }

    memset(data_, 0, elements_*sizeof(T));
    return 0;
  }

  T* get_data_ptr() { return data_; }
  
 private:
  std::vector<int> dimensions_;
  T* data_;
  long int elements_;
  
  int allocate_memory()
  {
    deallocate_memory();
    
    elements_ = 1;
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      elements_ *= dimensions_[i];
    } 
    
    try {
      data_ = new T[elements_];
    } catch (std::bad_alloc&) {
      std::cout << "NDArray<>::allocate memory failed" << std::endl;
      return -1;
    }
    
    clear();
    return 0;
  }
  
  int deallocate_memory() {
    if (data_) {
      delete [] data_;
      data_ = 0;
    }
    
    return 0;
  }

};


#endif //NDARRAY_H
