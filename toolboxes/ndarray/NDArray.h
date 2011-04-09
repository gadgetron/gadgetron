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

  virtual T* create(std::vector<unsigned int> dimensions) = 0;

  virtual T* create(std::vector<unsigned int> dimensions, T* data, 
		    bool delete_data_on_destruct = false) = 0;

  virtual int permute(std::vector<unsigned int>& dim_order, 
		      NDArray<T>* out = 0,
		      int shift_mode = 0) = 0;

  virtual int shift_dim(int shift, NDArray<T>* out = 0) {
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%dimensions_.size()));
    }
    return permute(order,out);
  }

  unsigned int get_number_of_dimensions() const {
    return dimensions_.size();
  }

  unsigned int get_size(unsigned int dimension) const {
    if (dimension >= dimensions_.size()) {
      return 1;
    } else {
      return dimensions_[dimension];
    }
  }

  std::vector<unsigned int> get_dimensions() const {
    return dimensions_;
  }

  unsigned long int get_number_of_elements() const {
    return elements_;
  }

  /**
     Remove all non-singleton dimensions;
   */
  void squeeze()
  {
    std::vector<unsigned int> new_dimensions;
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      if (dimensions_[i] != 1) {
	new_dimensions.push_back(dimensions_[i]);
      }
    }
    dimensions_ = new_dimensions;
  }

  int reshape(std::vector<unsigned int> dims) {
    unsigned long int new_elements = 1;
    for (unsigned int i = 0; i < dims.size(); i++) {
      new_elements *= dims[i];
    }

    if (new_elements != elements_) {
      std::cerr << "NDArray<T>::reshape : Number of elements cannot change during reshape" << std::endl;
      return -1;
    }
    
    dimensions_ = dims;
    return 0;

  }

  bool dimensions_equal(std::vector<unsigned int>& d) {
    return ((this->dimensions_.size() == d.size()) &&
	    std::equal(this->dimensions_.begin(), this->dimensions_.end(), d.begin()));
  }

  bool dimensions_equal(const NDArray<T>& a) {
    return ((this->dimensions_.size() == a.dimensions_.size()) &&
	    std::equal(this->dimensions_.begin(), this->dimensions_.end(), a.dimensions_.begin()));
  }

  T* get_data_ptr() { return data_; }
  
 protected:
  std::vector<unsigned int> dimensions_;
  T* data_;
  unsigned long int elements_;
  bool delete_data_on_destruct_;
  
};


#endif //NDARRAY_H
