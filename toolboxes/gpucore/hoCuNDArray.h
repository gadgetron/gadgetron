#pragma once
#include "hoNDArray.h"

template<class T> class hoCuNDArray: public hoNDArray<T>{

		hoCuNDArray() : hoNDArray<T>() {}

	  hoCuNDArray(std::vector<unsigned int> *dimensions) : hoNDArray<T>(dimensions) {}
	  hoCuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct=false) :
	  	hoNDArray<T> (dimensions, data, delete_data_on_destruct){};

	  hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions) : hoNDArray<T>(dimensions) {}
	  hoCuNDArray(boost::shared_ptr< std::vector<unsigned int> > dimensions, T* data, bool delete_data_on_destruct=false) :
		 hoNDArray<T> (dimensions, data, delete_data_on_destruct){};

};
