#ifndef HONDARRAY_H
#define HONDARRAY_H

#include "NDArray.h"

class ArrayIterator
{
 public:
  ArrayIterator(std::vector<unsigned int>dimensions, std::vector<unsigned int> order)
  {
    block_sizes_.push_back(1);
    for (unsigned int i = 0; i < order.size(); i++) {
      dimensions_.push_back(dimensions[i]);
      order_.push_back(order[i]);
      current_.push_back(0);
      if (i > 0) {
	block_sizes_.push_back(block_sizes_[i-1]*dimensions_[i-1]);
      }
    }
    current_idx_ = 0;
  }

  unsigned long int advance()
  {
    unsigned int order_index = 0;
    current_[order_[order_index]]++;
    while (current_[order_[order_index]] >= dimensions_[order_[order_index]]) {
      current_[order_[order_index]] = 0;
      order_index = (order_index+1)%dimensions_.size();
      current_[order_[order_index]]++;   
    }

    current_idx_ = 0;
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      current_idx_ += current_[i]*block_sizes_[i];
    }

    return current_idx_;
  }

  unsigned long int get_current_idx() 
  {
    return current_idx_;
  }

  std::vector<unsigned int>& get_current_sub()
  {
    return current_;
  }

 protected:
  std::vector<unsigned int> dimensions_;
  std::vector<unsigned int> order_;
  std::vector<unsigned int> current_;
  std::vector<unsigned long int> block_sizes_;

  unsigned long int current_idx_;
  
};

template <class T> class hoNDArray : public NDArray<T>
{
 public:
  hoNDArray () 
    : NDArray<T>::NDArray()
  {
    
  }

  ~hoNDArray() {
    if (this->delete_data_on_destruct_) {
      deallocate_memory();
    }
  }

  virtual T* create(std::vector<unsigned int>& dimensions) {
    this->dimensions_ = dimensions; 
    allocate_memory();
    return this->get_data_ptr();
  }

  virtual T* create(std::vector<unsigned int>& dimensions, T* data, 
		    bool delete_data_on_destruct = false) 
  {

    if (!data) {
      std::cerr << "hoNDArray::create : zero pointer provided" << std::endl;
      return 0;
    }

    this->dimensions_ = dimensions;
    this->data_ = data;
    this->delete_data_on_destruct_ = delete_data_on_destruct;
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
      this->elements_ *= this->dimensions_[i];
    }
    
    return this->get_data_ptr();
  }



  //Copy constructor
  hoNDArray(const hoNDArray<T>& a) {
    this->data_ = 0;
    this->dimensions_ = a.dimensions_;
    if (allocate_memory() == 0) {
      memcpy( this->data_, a.data_, this->elements_*sizeof(T) );
    }
  }
  
  
  virtual int permute(std::vector<unsigned int>& dim_order, NDArray<T>* out = 0)
  {
    hoNDArray<T>* out_int = 0;

    //Check ordering array
    if (dim_order.size() > this->dimensions_.size()) {
      std::cerr << "hoNDArray::permute - Invalid length of dimension ordering array" << std::endl;
      return -1;
    }

    std::vector<unsigned int> dim_count(this->dimensions_.size(),0);
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_order[i] >= this->dimensions_.size()) {
	std::cerr << "hoNDArray::permute - Invalid dimension order array" << std::endl;
	return -1;
      }
      dim_count[dim_order[i]]++;
    }

    //Create an internal array to store the dimensions
    std::vector<unsigned int> dim_order_int;

    //Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_count[dim_order[i]] != 1) {
	std::cerr << "hoNDArray::permute - Invalid dimension order array (duplicates)" << std::endl;
	return -1;
      }
      dim_order_int.push_back(dim_order[i]);
    }

    //Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < this->dimensions_.size()) {
      for (unsigned int i = 0; i < dim_count.size(); i++) {
	if (dim_count[i] == 0) {
	  dim_order_int.push_back(i);
	}
      }
    }

    if (out) {
      out_int = dynamic_cast< hoNDArray<T>* >(out);
      if (!out_int) {
	std::cerr << "hoNDArray::permute: failed to dynamic cast out array pointer" << std::endl;
	return -1;
      }
      for (unsigned int i = 0; i < dim_order_int.size(); i++) {
	if (this->dimensions_[dim_order_int[i]] != out_int->get_size(i)) {
	  std::cerr << "hoNDArray::permute: Dimensions of output array do not match the input array" << std::endl;
	  return -1;
	}
      }
    }


    T* o = 0;
    if (out_int) {
      o = out_int->get_data_ptr();
    } else {
      o = new T[this->elements_];
    }

    ArrayIterator it(this->dimensions_,dim_order_int);
    for (unsigned long int i = 0; i < this->elements_; i++) {
      o[i] = this->data_[it.get_current_idx()];
      it.advance();
    }

    if (!out_int) {
      std::vector<unsigned int> tmp_dims;
      for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
	tmp_dims.push_back(this->dimensions_[dim_order_int[i]]);
      }
      this->dimensions_ = tmp_dims;

      delete [] this->data_;
      this->data_ = o;
    }
   
    return 0;
  }

 protected:
  virtual int allocate_memory()
  {
    deallocate_memory();
    
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
      this->elements_ *= this->dimensions_[i];
    } 
    
    try {
      this->data_ = new T[this->elements_];
    } catch (std::bad_alloc&) {
      std::cout << "NDArray<>::allocate memory failed" << std::endl;
      return -1;
    }
    
    return 0;
  }
  
  virtual int deallocate_memory() {
    if (this->data_) {
      delete [] this->data_;
      this->data_ = 0;
    }
    
    return 0;
  }

};


#endif //HONDARRAY_H
