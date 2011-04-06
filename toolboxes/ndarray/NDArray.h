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
    , delete_data_on_destruct_(true)
  {
    
  }

  virtual ~NDArray() {  }

  virtual T* create(std::vector<unsigned int>& dimensions) = 0;

  virtual T* create(std::vector<unsigned int>& dimensions, T* data, 
		    bool delete_data_on_destruct = false) = 0;

  virtual int permute(std::vector<unsigned int>& dim_order, NDArray<T>* out = 0) = 0;

  virtual int shift_dim(int shift, NDArray<T>* out = 0) {
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%dimensions_.size()));
    }
    return permute(order,out);
  }

  unsigned int get_number_of_dimensions() {
    return dimensions_.size();
  }

  unsigned int get_size(unsigned int dimension) {
    if (dimension >= dimensions_.size()) {
      return 1;
    } else {
      return dimensions_[dimension];
    }
  }

  std::vector<unsigned int> get_dimensions() {
    return dimensions_;
  }

  unsigned long int get_number_of_elements() {
    return elements_;
  }

  T* get_data_ptr() { return data_; }
  
 protected:
  std::vector<unsigned int> dimensions_;
  T* data_;
  unsigned long int elements_;
  bool delete_data_on_destruct_;
  
};


#endif //NDARRAY_H
