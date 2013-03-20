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
#include "thrust/device_vector.h"
#include "GadgetronCuException.h"

namespace Gadgetron{
template <class T> class EXPORTGPUCORE cuNDArray;
template <class T> EXPORTGPUCORE void cuNDArray_permute(cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *order, int shift_mode);

template <class T> class EXPORTGPUCORE cuNDArray : public NDArray<T>
{
  friend void cuNDArray_permute<>(cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *order, int shift_mode);
    
 public:

  // Default constructor
  cuNDArray();
    
  // Copy constructor
  cuNDArray(const cuNDArray<T>& a);

  // Constructor from hoNDArray
  cuNDArray(hoNDArray<T> *a);
  cuNDArray(const hoNDArray<T>& a);

  cuNDArray(std::vector<unsigned int> *dimensions);
  cuNDArray(std::vector<unsigned int> *dimensions, int device_no);
  cuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions, int device_no);
  cuNDArray(boost::shared_ptr<std::vector<unsigned int>  > dimensions, T* data, bool delete_data_on_destruct = false);
  // Assignment operator
  cuNDArray& operator=(const cuNDArray<T>& rhs);
  
  virtual ~cuNDArray();

  virtual void create(std::vector<unsigned int> *dimensions);
  virtual void create(std::vector<unsigned int> *dimensions, int device_no);
  virtual void create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
  virtual void create(boost::shared_ptr<std::vector<unsigned int>  > dimensions){
  	this->create(dimensions.get());
  }
	virtual void create(boost::shared_ptr<std::vector<unsigned int>  > dimensions, int device_no){
		this->create(dimensions.get(),device_no);
	}
	virtual void create(boost::shared_ptr<std::vector<unsigned int>  > dimensions, T* data, bool delete_data_on_destruct = false){
		this->create(dimensions.get(), data, delete_data_on_destruct);
	}

  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions);
  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions, int device_no);  

  virtual boost::shared_ptr< hoNDArray<T> > to_host() const;
  virtual int to_host( hoNDArray<T>* ) const;

  virtual void permute(std::vector<unsigned int> *dim_order, NDArray<T> *out = 0, int shift_mode = 0);
  
  virtual void set_device(int device_no);
  inline int get_device() { return device_; }
  
  //Elementwise operations
  void abs();
  void sqrt();
  void clear();
  void fill(T );
  void reciprocal();
  void reciprocal_sqrt();


  thrust::device_ptr<T> get_device_ptr(){
	  return thrust::device_ptr<T>(this->data_);
  }

  thrust::device_ptr<T> begin(){
  	  return thrust::device_ptr<T>(this->data_);
  }
  thrust::device_ptr<T> end(){
  	  return thrust::device_ptr<T>(this->data_)+this->get_number_of_elements();
  }


 protected:
  
  int device_; 
  virtual void allocate_memory();
  virtual void deallocate_memory();

  
};

//Operators

 template<class T> void operator+= (cuNDArray<T> &, cuNDArray<T> &);
 template<class T> void operator+= (cuNDArray<T> &, T );

 template<class T> void operator*= (cuNDArray<T> &, cuNDArray<T> &);
 template<class T> void operator*= (cuNDArray<T> &, T );

 template<class T> void operator-= (cuNDArray<T> &, cuNDArray<T> &);
 template<class T> void operator-= (cuNDArray<T> & ,T );

 template<class T> void operator/= (cuNDArray<T> &, cuNDArray<T> &);
 template<class T> void operator/= (cuNDArray<T> &, T );

 template<class T> void operator+= (cuNDArray<complext<T> > &, cuNDArray<T> &);
template<class T> void operator+= (cuNDArray<complext<T> > &, T );

template<class T> void operator*= (cuNDArray<complext<T> > &, cuNDArray<T> &);
template<class T> void operator*= (cuNDArray<complext<T> > &, T );

template<class T> void operator-= (cuNDArray<complext<T> > &, cuNDArray<T> &);
template<class T> void operator-= (cuNDArray<complext<T> > & ,T );

template<class T> void operator/= (cuNDArray<complext<T> > &, cuNDArray<T> &);
template<class T> void operator/= (cuNDArray<complext<T> > &, T );


}


#endif //CUNDARRAY_H
