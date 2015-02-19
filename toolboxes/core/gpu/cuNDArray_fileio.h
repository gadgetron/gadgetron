#pragma once
#include "hoNDArray_fileio.h"

namespace Gadgetron{
template<class T> void write_nd_array(cuNDArray<T>* array, std::string s){
	write_nd_array(array->to_host().get(),s.c_str());
}


}
