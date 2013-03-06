#include "hoGTBLAS.h"
#include "gtArmadillo.h"

using namespace Gadgetron;
/**
 * @brief Returns the dot product of two arrays
 * @param[in] x
 * @param[in] y
 * @return dot product
 */
template<class T> T Gadgetron::dot(hoNDArray<T> *x,hoNDArray<T> *y){
	arma::Col<typename stdType<T>::type> xM = as_vector(x);
	arma::Col<typename stdType<T>::type> yM = as_vector(y);
	typename stdType<T>::type res =  arma::dot(xM,yM);

	return *(T*)(&res);
}

/**
 * @brief Returns the l2-norm of an array
 * @param[in] x
 * @return l2-norm
 */
template<class T> typename realType<T>::type Gadgetron::nrm2(hoNDArray<T> *x){
	typedef typename realType<T>::type real;
	arma::Col<typename stdType<T>::type> xM = as_vector(x);
	return real(arma::norm(xM,2));
}

/**
 * @brief Calculates Y = a*X+Y
 * @param[in] a A scalar
 * @param[in] x An array interpreted as a vector
 * @param[in,out] y An array interpreted as a vector
 */
template<class T> void Gadgetron::axpy(T a, hoNDArray<T>* x, hoNDArray<T>* y){

	typedef typename stdType<T>::type aT;
	arma::Col<aT> xM = as_vector(x);
	arma::Col<aT> yM = as_vector(y);
	aT a2 = *(aT*)(&a);
	yM = a2*xM+yM;

}
/**
 * @brief Calculates the index with the absolute minimum value
 * @param x
 * @return
 */
template<class T> int Gadgetron::amin(hoNDArray<T> *x){
	typename stdType<T>::type* ptr = (typename stdType<T>::type*)x->get_data_ptr();
	typename realType<T>::type cur_min = std::abs(ptr[0]);
	int index = 0;
	for(int i=1; i< x->get_number_of_elements(); i++)
	{
		typename realType<T>::type val = std::abs(ptr[i]);
		if (val < cur_min){
			cur_min = val;
			index = i;
		}
	}
	return index;
}

/**
 * @brief Calculates the index with the absolute maximum value
 * @param x
 * @return
 */
template<class T> int Gadgetron::amax(hoNDArray<T> *x){

	typename stdType<T>::type* ptr = (typename stdType<T>::type*)x->get_data_ptr();
	typename realType<T>::type cur_max = std::abs(ptr[0]);
	int index = 0;
	for(int i=1; i< x->get_number_of_elements(); i++)
	{
		typename realType<T>::type val = std::abs(ptr[i]);
		if (val > cur_max){
			cur_max = val;
			index = i;
		}
	}
	return index;
}

//TODO: Make sure this approach doesn't actually copy our x....
// 01-03-2013: Runtime seems to be the same as the simple loop approach, so it's probably safe
template<class T> typename realType<T>::type Gadgetron::asum(hoNDArray<T> *x){
	arma::Col<typename stdType<T>::type> xM = as_vector(x);
	return arma::accu(arma::abs(xM));

}


template double Gadgetron::dot(hoNDArray<double>*,hoNDArray<double>*);
template double Gadgetron::nrm2(hoNDArray<double> *x);
template void Gadgetron::axpy(double a, hoNDArray<double>* x, hoNDArray<double>* y);
template int Gadgetron::amin(hoNDArray<double>* x);
template int Gadgetron::amax(hoNDArray<double>* x);
template double Gadgetron::asum(hoNDArray<double>* x);

template float Gadgetron::dot(hoNDArray<float>*,hoNDArray<float>*);
template float Gadgetron::nrm2(hoNDArray<float> *x);
template void Gadgetron::axpy(float a, hoNDArray<float>* x, hoNDArray<float>* y);
template int Gadgetron::amin(hoNDArray<float>* x);
template int Gadgetron::amax(hoNDArray<float>* x);
template float Gadgetron::asum(hoNDArray<float>* x);

template double_complext Gadgetron::dot(hoNDArray<double_complext>*,hoNDArray<double_complext>*);
template double Gadgetron::nrm2(hoNDArray<double_complext> *x);
template void Gadgetron::axpy(double_complext a, hoNDArray<double_complext>* x, hoNDArray<double_complext>* y);
template int Gadgetron::amin(hoNDArray<double_complext>* x);
template int Gadgetron::amax(hoNDArray<double_complext>* x);
template double Gadgetron::asum(hoNDArray<double_complext>* x);

template float_complext Gadgetron::dot(hoNDArray<float_complext>*,hoNDArray<float_complext>*);
template float Gadgetron::nrm2(hoNDArray<float_complext> *x);
template void Gadgetron::axpy(float_complext a, hoNDArray<float_complext>* x, hoNDArray<float_complext>* y);
template int Gadgetron::amin(hoNDArray<float_complext>* x);
template int Gadgetron::amax(hoNDArray<float_complext>* x);
template float Gadgetron::asum(hoNDArray<float_complext>* x);
