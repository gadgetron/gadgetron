#ifndef CUNDARRAY_H
#define CUNDARRAY_H
#pragma once

#include "gpucore_export.h"
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "NDArray.h"
#include "hoNDArray.h"

#include <boost/shared_ptr.hpp>
#include "complext.h"

template <class T> class EXPORTGPUCORE cuNDArray;
template <class T> EXPORTGPUCORE int cuNDArray_permute(cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *order, int shift_mode);

template <class T> class EXPORTGPUCORE cuNDArray : public NDArray<T>
{
  friend int cuNDArray_permute<>(cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *order, int shift_mode);
    
 public:

  // Default constructor
  cuNDArray();
    
  // Copy constructor
  cuNDArray(const cuNDArray<T>& a);

  // Constructor from hoNDArray
  cuNDArray(hoNDArray<T> *a);

  cuNDArray(std::vector<unsigned int> *dimensions);
  cuNDArray(std::vector<unsigned int> *dimensions, int device_no);
  cuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions, int device_no);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions, T* data, bool delete_data_on_destruct = false);
  // Assignment operator
  cuNDArray& operator=(const cuNDArray<T>& rhs);
  
  virtual ~cuNDArray();

  virtual T* create(std::vector<unsigned int> *dimensions);
  virtual T* create(std::vector<unsigned int> *dimensions, int device_no);
  virtual T* create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);

  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions);
  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions, int device_no);  

  virtual boost::shared_ptr< hoNDArray<T> > to_host() const;

  virtual int permute(std::vector<unsigned int> *dim_order, NDArray<T> *out = 0, int shift_mode = 0);
  
  virtual int set_device(int device_no);
  inline int get_device() { return device_; }
  
 protected:
  
  int device_; 
  virtual int allocate_memory();
  virtual int deallocate_memory();
  
};

#endif //CUNDARRAY_H
