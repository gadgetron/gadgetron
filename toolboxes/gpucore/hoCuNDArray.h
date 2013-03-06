#pragma once
#include "hoNDArray.h"

namespace Gadgetron{
template<class T> class hoCuNDArray: public hoNDArray<T>{
public:
		hoCuNDArray() : hoNDArray<T>() {}

	  hoCuNDArray(std::vector<unsigned int> *dimensions) : hoNDArray<T>(dimensions) {}
	  hoCuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct=false) :
	  	hoNDArray<T> (dimensions, data, delete_data_on_destruct){};

	  hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions) : hoNDArray<T>(dimensions) {}
	  hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions, T* data, bool delete_data_on_destruct=false) :
		 hoNDArray<T> (dimensions, data, delete_data_on_destruct){};

	  // Copy constructor
		hoCuNDArray(const hoNDArray<T>& a)
		{
	      this->data_ = 0;
	      this->dimensions_ = a.dimensions_;
	      this->allocate_memory();
	      memcpy( this->data_, a.data_, this->elements_*sizeof(T) );

		}

	    // Assignment operator
		hoCuNDArray& operator=(const hoNDArray<T>& rhs) {
			// Are the dimensions the same? Then we can just memcpy
			if (this->dimensions_equal(&rhs)) {
				memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
			} else {
				this->deallocate_memory();
				this->data_ = 0;
				this->dimensions_ = rhs.dimensions_;
				this->allocate_memory();
				memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );

			}
			return *this;
		}

};
}
