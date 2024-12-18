#pragma once
#include "complext.h"
#include "hoCuNDArray.h"
#include "cuNDArray.h"

namespace Gadgetron{

template<class T> void solver_non_negativity_filter(cuNDArray<T>* x , cuNDArray<T>* g);


template<class T> void updateF(cuNDArray<T>& data, typename realType<T>::Type alpha ,typename realType<T>::Type sigma);

template<class T> void updateF(hoCuNDArray<T>& data, typename realType<T>::Type alpha ,typename realType<T>::Type sigma){
	cuNDArray<T> cudata(data);
	updateF<T>(cudata,alpha,sigma);
	data = cudata;

}

template<class T> void updateFgroup(std::vector<cuNDArray<T> >& datas, typename realType<T>::Type alpha ,typename realType<T>::Type sigma);
template<class T> void updateFgroup(std::vector<hoCuNDArray<T> >& datas, typename realType<T>::Type alpha ,typename realType<T>::Type sigma){
	std::vector<cuNDArray<T> > cudatas(datas.size());
	for (size_t i =0; i < datas.size(); i++)
		cudatas[i] = cuNDArray<T>(datas[i]);

	updateFgroup<T>(cudatas,alpha,sigma);

	for (size_t i =0; i < datas.size(); i++)
		datas[i] = cudatas[i];

}

}
