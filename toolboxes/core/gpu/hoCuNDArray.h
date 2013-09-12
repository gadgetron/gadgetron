/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container) for cpu->gpu->cpu (hoCu) solvers.

    The existence of this class is mainly due to providing unique array type for the hoCu based math in
    hoCuNDArray_operators.h, hoCuNDArray_elemwise.h, and hoCuNDArray_blas.h.
    Unfortunately C++ does not let a derived class inherit its base class's constructors, which consequently need redefinition.
 */

#pragma once

#include "hoNDArray.h"
#include "check_CUDA.h"

namespace Gadgetron{

template<class T> class hoCuNDArray: public hoNDArray<T>
{
public:

	hoCuNDArray() : hoNDArray<T>::hoNDArray() {}

	hoCuNDArray(std::vector<unsigned int> *dimensions) : hoNDArray<T>::hoNDArray() {
		this->create(dimensions);
	}

	hoCuNDArray(std::vector<unsigned int> *dimensions,
			T* data, bool delete_data_on_destruct = false) : hoNDArray<T>::hoNDArray() {
		this->create(dimensions,data,delete_data_on_destruct);
	}

	hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions) : hoNDArray<T>::hoNDArray() {
		this->create(dimensions.get());
	}

	hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions,
			T* data, bool delete_data_on_destruct = false) : hoNDArray<T>::hoNDArray() {
		this->create(dimensions.get(),data,delete_data_on_destruct);
	}

	// Copy constructor
	hoCuNDArray(const hoNDArray<T>& a): hoNDArray<T>(){
		this->data_ = 0;
		this->dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(new std::vector<unsigned int>(*a.get_dimensions()));
		this->allocate_memory();
		memcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T));
	}

	// Copy constructor
	hoCuNDArray(const hoCuNDArray<T>& a): hoNDArray<T>(){
		this->data_ = 0;
		this->dimensions_ = boost::shared_ptr< std::vector<unsigned int> >(new std::vector<unsigned int>(*a.get_dimensions()));
		this->allocate_memory();
		memcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T));
	}
	virtual ~hoCuNDArray() {
		if (this->delete_data_on_destruct_) {
			this->deallocate_memory();
		}
	}

protected:

	virtual void allocate_memory()
	{

		this->deallocate_memory();
		this->elements_ = 1;
		if (this->dimensions_->empty())
			throw std::runtime_error("hoCuNDArray::allocate_memory() : dimensions is empty.");
		for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
			this->elements_ *= (*this->dimensions_)[i];
		}

		size_t size = this->elements_ * sizeof(T);

		//this->data_= (T*)std::malloc(size);
		CUDA_CALL(cudaMallocHost((void**)&this->data_,size));
	}

	virtual void deallocate_memory()
	{
		if (this->data_) {
			CUDA_CALL(cudaFreeHost(this->data_));
			//free(this->data_);
			this->data_ = 0;
		}
	}



};

}
