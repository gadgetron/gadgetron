#include "hoNDArray_reductions.h"
#include "hoArmadillo.h"

namespace Gadgetron{

template<class REAL> REAL max(hoNDArray<REAL>* data){
	return as_arma_col(data).max();
}
template<class REAL> REAL min(hoNDArray<REAL>* data){
	return as_arma_col(data).min();
}


template<class T> T mean(hoNDArray<T>* data){
	return (typename stdType<T>::Type) arma::mean(as_arma_col(data));
}


template<class T> T sum(hoNDArray<T>* data){
	return (typename stdType<T>::Type) arma::sum(as_arma_col(data));
}


template EXPORTCPUCOREMATH float max(hoNDArray<float>*);
template EXPORTCPUCOREMATH float min(hoNDArray<float>*);
template EXPORTCPUCOREMATH float mean(hoNDArray<float>*);
template EXPORTCPUCOREMATH float sum(hoNDArray<float>*);

template EXPORTCPUCOREMATH double max(hoNDArray<double>*);
template EXPORTCPUCOREMATH double min(hoNDArray<double>*);
template EXPORTCPUCOREMATH double mean(hoNDArray<double>*);
template EXPORTCPUCOREMATH double sum(hoNDArray<double>*);


template EXPORTCPUCOREMATH complext<double> mean(hoNDArray<complext<double> >*);
template EXPORTCPUCOREMATH complext<double> sum(hoNDArray<complext<double> >*);

template EXPORTCPUCOREMATH complext<float> mean(hoNDArray<complext<float> >*);
template EXPORTCPUCOREMATH complext<float> sum(hoNDArray<complext<float> >*);
}

