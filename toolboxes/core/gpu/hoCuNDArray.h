/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container) for cpu->gpu->cpu (hoCu) solvers.

    The existence of this class is mainly due to providing unique array type for the hoCu based math in
    hoCuNDArray_operators.h, hoCuNDArray_elemwise.h, and hoCuNDArray_blas.h.
    Unfortunately C++ does not let a derived class inherit its base class's constructors, which consequently need redefinition.
*/

#pragma once

#include "cuNDArray.h"
#include "hoNDArray.h"
#include "check_CUDA.h"

namespace Gadgetron{

  template<class T> class hoCuNDArray: public hoNDArray<T>
  {
  public:

    hoCuNDArray() : Gadgetron::hoNDArray<T>::hoNDArray() {}


#if __cplusplus > 199711L
    hoCuNDArray(hoCuNDArray<T>&& other) : Gadgetron::hoNDArray<T>::hoNDArray(){
        this->data_ = other.data_;
        this->dimensions_ = other.dimensions_;
        this->elements_ = other.elements_;
        this->offsetFactors_ = std::move(other.offsetFactors_);
        other.data_ = nullptr;
    }
#endif

    hoCuNDArray(std::vector<size_t> *dimensions) : Gadgetron::hoNDArray<T>::hoNDArray() {
      this->create(dimensions);
    }

    hoCuNDArray(std::vector<size_t> &dimensions) : Gadgetron::hoNDArray<T>::hoNDArray() {
      this->create(dimensions);
    }
  
    hoCuNDArray(std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false) : Gadgetron::hoNDArray<T>::hoNDArray() {
      this->create(dimensions, data, delete_data_on_destruct);
    }
  
    // Copy constructors
    hoCuNDArray(const hoNDArray<T> &a): Gadgetron::hoNDArray<T>::hoNDArray(){
      this->create(a.get_dimensions());
      memcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T));
    }

    hoCuNDArray(const hoNDArray<T> *a): Gadgetron::hoNDArray<T>::hoNDArray(){
      if(!a) throw std::runtime_error("hoCuNDArray::hoCuNDArray(): 0x0 pointer provided.");
      this->create(a->get_dimensions());
      memcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T));
    }

    hoCuNDArray(const hoCuNDArray<T> &a): Gadgetron::hoNDArray<T>::hoNDArray(){
      this->create(a.get_dimensions());
      memcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T));
    }

    hoCuNDArray(const hoCuNDArray<T> *a): Gadgetron::hoNDArray<T>::hoNDArray(){
      if(!a) throw std::runtime_error("hoCuNDArray::hoCuNDArray(): 0x0 pointer provided.");
      this->create(a->get_dimensions());
      memcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T));
    }

    virtual ~hoCuNDArray() {
      if (this->delete_data_on_destruct_) {
        this->deallocate_memory();
      }
    }

    T& at( size_t idx ){
      if( idx >= this->get_number_of_elements() ){
        throw std::runtime_error("hoCuNDArray::at(): index out of range.");
      }
      return this->data_[idx];
    }
  
    T& operator[]( size_t idx ){
      if( idx >= this->get_number_of_elements() ){
        throw std::runtime_error("hoCuNDArray::operator[]: index out of range.");
      }
      return this->data_[idx];
    }

    hoCuNDArray<T>& operator=(const hoCuNDArray<T>& rhs)
    {
        if ( &rhs == this ) return *this;

        if ( rhs.get_number_of_elements() == 0 ){
            this->clear();
            return *this;
        }

        // Are the dimensions the same? Then we can just memcpy
        if (this->dimensions_equal(rhs)){
            memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
        }
        else{
            deallocate_memory();
            this->data_ = 0;
            this->dimensions_ = rhs.dimensions_;
            this->offsetFactors_ = rhs.offsetFactors_;
            this->allocate_memory();
            memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );
        }
        return *this;
    }


    hoCuNDArray<T>& operator=(const cuNDArray<T> & rhs)
    {

        if ( rhs.get_number_of_elements() == 0 ){
            this->clear();
            return *this;
        }

        // Are the dimensions the same? Then we can just memcpy
        if (!this->dimensions_equal(rhs)){
            this->create(rhs.get_dimensions());
        }

          cudaMemcpy(this->data_, rhs.get_data_ptr(), this->elements_*sizeof(T),cudaMemcpyDeviceToHost);
        return *this;
    }
#if __cplusplus > 199711L
    hoCuNDArray<T>& operator=(hoCuNDArray<T>&& rhs)
    {
        if ( &rhs == this ) return *this;

        this->clear();
        this->dimensions_ = rhs.dimensions_;
        this->offsetFactors_ = rhs.offsetFactors_;
        this->elements_ = rhs.elements_;
        this->data_ = rhs.data_;
        rhs.data_ = nullptr;
        return *this;
    }
#endif
  protected:

    virtual void allocate_memory()
    {
      this->deallocate_memory();
      this->elements_ = 1;
      if (this->dimensions_.empty())
        throw std::runtime_error("hoCuNDArray::allocate_memory() : dimensions is empty.");
      for (size_t i = 0; i < this->dimensions_.size(); i++) {
        this->elements_ *= this->dimensions_[i];
      }

      size_t size = this->elements_ * sizeof(T);
      CUDA_CALL(cudaMallocHost((void**)&this->data_,size));
    }

    virtual void deallocate_memory()
    {
      if (this->data_) {
        CUDA_CALL(cudaFreeHost(this->data_));
        this->data_ = 0;
      }
    }
  };
}
