/** \file cuNDArray.h
    \brief GPU-based N-dimensional array (data container)
*/

#ifndef CUNDARRAY_H
#define CUNDARRAY_H
#pragma once

#include "NDArray.h"
#include "hoNDArray.h"
#include "complext.h"
#include "GadgetronCuException.h"
#include "check_CUDA.h"

#include <boost/shared_ptr.hpp>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>

namespace Gadgetron{

  template <class T> class cuNDArray : public NDArray<T>
  {

  public:

    // Constructors
    //

    cuNDArray() : NDArray<T>::NDArray() 
    { 
      cudaGetDevice(&this->device_); 
    }

    cuNDArray(const cuNDArray<T> &a) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      this->data_ = 0;
      this->dimensions_ = a.get_dimensions();
      allocate_memory();
      if (a.device_ == this->device_) {
        CUDA_CALL(cudaMemcpy(this->data_, a.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice));
      } else {
        //This memory is on a different device, we must move it.
        cudaSetDevice(a.device_);
        boost::shared_ptr< hoNDArray<T> > tmp = a.to_host();
        cudaSetDevice(this->device_);
        cudaError_t err = cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice);
        if (err !=cudaSuccess) {
          deallocate_memory();
          this->data_ = 0;
          this->dimensions_->clear();
          throw cuda_error(err);
        }
      }
    }

    cuNDArray(const cuNDArray<T> *a) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      this->data_ = 0;
      this->dimensions_ = a->get_dimensions();
      allocate_memory();
      if (a->device_ == this->device_) {
        CUDA_CALL(cudaMemcpy(this->data_, a->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice));
      } else {
        //This memory is on a different device, we must move it.
        cudaSetDevice(a->device_);
        boost::shared_ptr< hoNDArray<T> > tmp = a->to_host();
        cudaSetDevice(this->device_);
        cudaError_t err = cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice);
        if (err !=cudaSuccess) {
          deallocate_memory();
          this->data_ = 0;
          this->dimensions_->clear();
          throw cuda_error(err);
        }
      }
    }

    cuNDArray(const hoNDArray<T> &a) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      this->dimensions_ = a.get_dimensions();
      allocate_memory();
      if (cudaMemcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
        deallocate_memory();
        this->data_ = 0;
        this->dimensions_->clear();
      }
    }

    cuNDArray(hoNDArray<T> *a) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      this->dimensions_ = a->get_dimensions();
      allocate_memory();
      if (cudaMemcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
        deallocate_memory();
        this->data_ = 0;
        this->dimensions_->clear();
      }
    }

    cuNDArray(std::vector<size_t> *dimensions) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      create(dimensions);
    }

    cuNDArray(std::vector<size_t> *dimensions, int device_no) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      create(dimensions,device_no);
    }

    cuNDArray(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false) : NDArray<T>::NDArray()
    {
      cudaGetDevice(&this->device_);
      create(dimensions,data,delete_data_on_destruct);
    }

    cuNDArray(std::vector<size_t> &dimensions) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      create(dimensions);
    }

    cuNDArray(std::vector<size_t> &dimensions, int device_no) : NDArray<T>::NDArray() 
    {
      cudaGetDevice(&this->device_);
      create(dimensions,device_no);
    }

    cuNDArray(std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false) : NDArray<T>::NDArray()
    {
      cudaGetDevice(&this->device_);
      create(dimensions,data,delete_data_on_destruct);
    }

    cuNDArray(boost::shared_ptr<std::vector<size_t> > dimensions) : NDArray<T>::NDArray()
    {
      cudaGetDevice(&this->device_);
      create(dimensions.get());
    }

    cuNDArray(boost::shared_ptr<std::vector<size_t> > dimensions, int device_no) : NDArray<T>::NDArray()
    {
      cudaGetDevice(&this->device_);
      create(dimensions.get(),device_no);
    }

    cuNDArray(boost::shared_ptr<std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false) : NDArray<T>::NDArray()
    {
      cudaGetDevice(&this->device_);
      create(dimensions.get(),data,delete_data_on_destruct);
    }

    cuNDArray(size_t len)
    {
      std::vector<size_t> dim(1);
      dim[0] = len;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy)
    {
      std::vector<size_t> dim(2);
      dim[0] = sx;
      dim[1] = sy;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz)
    {
      std::vector<size_t> dim(3);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz, size_t st)
    {
      std::vector<size_t> dim(4);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp)
    {
      std::vector<size_t> dim(5);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq)
    {
      std::vector<size_t> dim(6);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr)
    {
      std::vector<size_t> dim(7);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      dim[6] = sr;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    cuNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss)
    {
      std::vector<size_t> dim(8);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      dim[6] = sr;
      dim[7] = ss;
      cudaGetDevice(&this->device_);
      create(dim);
    }

    // Destructor
    virtual ~cuNDArray()
    { 
      if (this->delete_data_on_destruct_) 
        deallocate_memory();  
    }

    // Assignment operator
    cuNDArray<T>& operator=(const cuNDArray<T>& rhs)
    {
      int cur_device; 
      CUDA_CALL(cudaGetDevice(&cur_device));
      bool dimensions_match = this->dimensions_equal(&rhs);
      if (dimensions_match && (rhs.device_ == cur_device) && (cur_device == this->device_)) {
        CUDA_CALL(cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice));
      }
      else {
        CUDA_CALL(cudaSetDevice(this->device_));
        if( !dimensions_match ){
          deallocate_memory();
          this->elements_ = rhs.elements_;
          this->dimensions_ = rhs.get_dimensions();
          allocate_memory();
        }
        if (this->device_ == rhs.device_) {
          if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=cudaSuccess) {	    
            cudaSetDevice(cur_device);
            throw cuda_error("cuNDArray::operator=: failed to copy data (2)");
          }
        } else {
          if( cudaSetDevice(rhs.device_) != cudaSuccess) {
            cudaSetDevice(cur_device);
            throw cuda_error("cuNDArray::operator=: unable to set device no (2)");
          }
          boost::shared_ptr< hoNDArray<T> > tmp = rhs.to_host();
          if( cudaSetDevice(this->device_) != cudaSuccess) {
            cudaSetDevice(cur_device);
            throw cuda_error("cuNDArray::operator=: unable to set device no (3)");
          }
          if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
            cudaSetDevice(cur_device);
            throw cuda_error("cuNDArray::operator=: failed to copy data (3)");
          }
        }
        if( cudaSetDevice(cur_device) != cudaSuccess) {
          throw cuda_error("cuNDArray::operator=: unable to restore to current device");
        }
      }
      return *this;
    }

    cuNDArray<T>& operator=(const hoNDArray<T>& rhs)
    {
      int cur_device; 
      CUDA_CALL(cudaGetDevice(&cur_device));
      bool dimensions_match = this->dimensions_equal(&rhs);
      if (dimensions_match && (cur_device == this->device_)) {
        CUDA_CALL(cudaMemcpy(this->get_data_ptr(), rhs.get_data_ptr(), this->get_number_of_elements()*sizeof(T), cudaMemcpyHostToDevice));
      }
      else {
        CUDA_CALL(cudaSetDevice(this->device_));
        if( !dimensions_match ){
          deallocate_memory();
          this->elements_ = rhs.get_number_of_elements();
          this->dimensions_ = rhs.get_dimensions();
          allocate_memory();
        }
        if (cudaMemcpy(this->get_data_ptr(), rhs.get_data_ptr(), this->get_number_of_elements()*sizeof(T),
                       cudaMemcpyHostToDevice) !=cudaSuccess) {
          cudaSetDevice(cur_device);
          throw cuda_error("cuNDArray::operator=: failed to copy data (1)");
        }
        if( cudaSetDevice(cur_device) != cudaSuccess) {
          throw cuda_error("cuNDArray::operator=: unable to restore to current device");
        }
      }
      return *this;
    }

    virtual void create(std::vector<size_t> *dimensions)
    {
      if ( this->dimensions_equal(dimensions) )
        {
          return;
        }

      return NDArray<T>::create(dimensions);
    }

    virtual void create(std::vector<size_t> *dimensions, int device_no)
    {
      if (device_no < 0){
        throw cuda_error("cuNDArray::create: illegal device no");
      }

      if ( this->dimensions_equal(dimensions) && this->device_==device_no )
        {
          return;
        }

      this->device_ = device_no; 
      NDArray<T>::create(dimensions);
    }

    virtual void create(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false)
    {
      if (!data) {
        throw std::runtime_error("cuNDArray::create: 0x0 pointer provided");
      }

      int tmp_device; 
      if( cudaGetDevice(&tmp_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::create: Unable to query for device");
      }

      cudaDeviceProp deviceProp; 
      if( cudaGetDeviceProperties( &deviceProp, tmp_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::create: Unable to query device properties");
      }

      if (deviceProp.unifiedAddressing) {
        cudaPointerAttributes attrib;
        if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
          CHECK_FOR_CUDA_ERROR();
          throw cuda_error("cuNDArray::create: Unable to determine attributes of pointer");
        }
        this->device_ = attrib.device;
      } else {
        this->device_ = tmp_device;
      }

      NDArray<T>::create(dimensions, data, delete_data_on_destruct);
    }

    virtual void create(std::vector<size_t> &dimensions)
    {
      if ( this->dimensions_equal(&dimensions) )
        {
          return;
        }

      return NDArray<T>::create(dimensions);
    }

    virtual void create(std::vector<size_t> &dimensions, int device_no)
    {
      if (device_no < 0){
        throw cuda_error("cuNDArray::create: illegal device no");
      }

      if ( this->dimensions_equal(&dimensions) && this->device_==device_no )
        {
          return;
        }

      this->device_ = device_no; 
      NDArray<T>::create(dimensions);
    }

    virtual void create(std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false)
    {
      if (!data) {
        throw std::runtime_error("cuNDArray::create: 0x0 pointer provided");
      }

      int tmp_device; 
      if( cudaGetDevice(&tmp_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::create: Unable to query for device");
      }

      cudaDeviceProp deviceProp;
      if( cudaGetDeviceProperties( &deviceProp, tmp_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::create: Unable to query device properties");
      }

      if (deviceProp.unifiedAddressing) {
        cudaPointerAttributes attrib;
        if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
          CHECK_FOR_CUDA_ERROR();
          throw cuda_error("cuNDArray::create: Unable to determine attributes of pointer");
        }
        this->device_ = attrib.device;
      } else {
        this->device_ = tmp_device;
      }

      NDArray<T>::create(dimensions, data, delete_data_on_destruct);
    }

    virtual void create(boost::shared_ptr<std::vector<size_t> > dimensions){
      this->create(dimensions.get());
    }

    virtual void create(boost::shared_ptr<std::vector<size_t> > dimensions, int device_no){
      this->create(dimensions.get(),device_no);
    }

    virtual void create(boost::shared_ptr<std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false){
      this->create(dimensions.get(), data, delete_data_on_destruct);
    }

    virtual void create(size_t len)
    {
      std::vector<size_t> dim(1);
      dim[0] = len;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy)
    {
      std::vector<size_t> dim(2);
      dim[0] = sx;
      dim[1] = sy;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz)
    {
      std::vector<size_t> dim(3);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz, size_t st)
    {
      std::vector<size_t> dim(4);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp)
    {
      std::vector<size_t> dim(5);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq)
    {
      std::vector<size_t> dim(6);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr)
    {
      std::vector<size_t> dim(7);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      dim[6] = sr;
      this->create(dim);
    }

    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss)
    {
      std::vector<size_t> dim(8);
      dim[0] = sx;
      dim[1] = sy;
      dim[2] = sz;
      dim[3] = st;
      dim[4] = sp;
      dim[5] = sq;
      dim[6] = sr;
      dim[7] = ss;
      this->create(dim);
    }

    virtual boost::shared_ptr< hoNDArray<T> > to_host() const
    {
      boost::shared_ptr< hoNDArray<T> > ret(new hoNDArray<T>(this->dimensions_.get()));
      if (cudaMemcpy(ret->get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw cuda_error("cuNDArray::to_host(): failed to copy memory from device");
      }

      return ret;
    }

    virtual void to_host( hoNDArray<T> *out ) const 
    {
      if( !out ){
        throw std::runtime_error("cuNDArray::to_host(): illegal array passed.");
      }

      if( out->get_number_of_elements() != this->get_number_of_elements() ){	
        out->create( this->get_dimensions().get());
      }

      if( cudaMemcpy( out->get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw cuda_error("cuNDArray::to_host(): failed to copy memory from device");
      }
    }

    virtual void set_device(int device)
    {
      if( device_ == device )
        return;

      int cur_device;
      if( cudaGetDevice(&cur_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::set_device: unable to get device no");
      }

      if( cur_device != device_ && cudaSetDevice(device_) != cudaSuccess) {
        throw cuda_error("cuNDArray::set_device: unable to set device no");
      }

      boost::shared_ptr< hoNDArray<T> > tmp = to_host();
      deallocate_memory();
      if( cudaSetDevice(device) != cudaSuccess) {
        cudaSetDevice(cur_device);
        throw cuda_error("cuNDArray::set_device: unable to set device no (2)");
      }

      device_ = device;
      allocate_memory();
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaSetDevice(cur_device);
        throw cuda_error("cuNDArray::set_device: failed to copy data");
      }

      if( cudaSetDevice(cur_device) != cudaSuccess) {
        throw cuda_error("cuNDArray::set_device: unable to restore device to current device");
      }
    }

    inline int get_device() { return device_; }

    thrust::device_ptr<T> get_device_ptr(){
      return thrust::device_ptr<T>(this->data_);
    }

    thrust::device_ptr<T> begin(){
      return thrust::device_ptr<T>(this->data_);
    }

    thrust::device_ptr<T> end(){
      return thrust::device_ptr<T>(this->data_)+this->get_number_of_elements();
    }

    T at( size_t idx ){
      if( idx >= this->get_number_of_elements() ){
        throw std::runtime_error("cuNDArray::at(): index out of range.");
      }
      T res;
      CUDA_CALL(cudaMemcpy(&res, &this->get_data_ptr()[idx], sizeof(T), cudaMemcpyDeviceToHost));
      return res;
    }

    T operator[]( size_t idx ){
      if( idx >= this->get_number_of_elements() ){
        throw std::runtime_error("cuNDArray::operator[]: index out of range.");
      }
      T res;
      CUDA_CALL(cudaMemcpy(&res, &this->get_data_ptr()[idx], sizeof(T), cudaMemcpyDeviceToHost));
      return res;
    }

  protected:

    int device_; 

    virtual void allocate_memory()
    {
      deallocate_memory();

      this->elements_ = 1;
      if (this->dimensions_->empty())
        throw std::runtime_error("cuNDArray::allocate_memory() : dimensions is empty.");
      for (size_t i = 0; i < this->dimensions_->size(); i++) {
        this->elements_ *= (*this->dimensions_)[i];
      } 

      size_t size = this->elements_ * sizeof(T);

      int device_no_old;
      if (cudaGetDevice(&device_no_old) != cudaSuccess) {
        throw cuda_error("cuNDArray::allocate_memory: unable to get device no");
      }

      if (device_ != device_no_old) {
        if (cudaSetDevice(device_) != cudaSuccess) {
          throw cuda_error("cuNDArray::allocate_memory: unable to set device no");
        }
      }

      if (cudaMalloc((void**) &this->data_,size) != cudaSuccess) {
        size_t free = 0, total = 0;
        cudaMemGetInfo(&free, &total);
        std::stringstream err("cuNDArray::allocate_memory() : Error allocating CUDA memory");
        err << "CUDA Memory: " << free << " (" << total << ")";

        err << "   memory requested: " << size << "( ";
        for (size_t i = 0; i < this->dimensions_->size(); i++) {
          std::cerr << (*this->dimensions_)[i] << " ";
        }
        err << ")";
        this->data_ = 0;
        throw std::runtime_error(err.str());
      }

      if (device_ != device_no_old) {
        if (cudaSetDevice(device_no_old) != cudaSuccess) {
          throw cuda_error("cuNDArray::allocate_memory: unable to restore device no");
        }
      }
    }

    virtual void deallocate_memory()
    {
      if (this->data_) {

        int device_no_old;
        CUDA_CALL(cudaGetDevice(&device_no_old));
        if (device_ != device_no_old) {
          CUDA_CALL(cudaSetDevice(device_));
        }

        CUDA_CALL(cudaFree(this->data_));
        if (device_ != device_no_old) {
          CUDA_CALL(cudaSetDevice(device_no_old));
        }
        this->data_ = 0;
      }
    }
  };
}

#endif //CUNDARRAY_H
