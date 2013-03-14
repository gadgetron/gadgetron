#pragma once

#include "gtArmadillo.h"

namespace Gadgetron{


/**
 * @brief Calculates the elementwise absolute value of the array
 * @param[in] data Input data
 * @return A new array containing the elementwise absolute value of data
 */
template<class T>
boost::shared_ptr<hoNDArray<typename realType<T>::type> > abs(hoNDArray<T> *data){
	typedef typename realType<T>::type REAL;
	typedef typename stdType<T>::type STD;
	boost::shared_ptr<hoNDArray<REAL> > res(new hoNDArray<REAL>(data->get_dimensions()));

	arma::Col<REAL> res_vec = as_vector(res.get());
	arma::Col<STD> data_vec = as_vector(data);

	res_vec = arma::abs(data_vec);
	return res;
}


/**
 * @brief Clamps all values in the array to the minimum and maximum values specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] min minimum value
 * @param[in] max maximum value
 */
template<class  T>
void clamp(hoNDArray<T> *in_out, T min, T max){
	T* x = in_out->get_data_ptr();
	for (int i = 0; i < in_out->get_number_of_elements(); i++){
		 if (x[i] < min) x[i]= min;
		 else if (x[i] > max) x[i] =  max;

	}
}


/**
 * @brief Clamps all values in the array to the minimum value specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] min minimum value
 */
template<class  T>
void clamp_min(hoNDArray<T> *in_out, T min){
	T* x = in_out->get_data_ptr();
	for (int i = 0; i < in_out->get_number_of_elements(); i++){
		 if (x[i] < min) x[i]= min;

	}
}

/**
 * @brief Clamps all values in the array to the maximum value specified.
 * @param[in,out] in_out Array which to clamp
 * @param[in] max minimum value
 */
template<class  T>
void clamp_max(hoNDArray<T> *in_out, T max){
	T* x = in_out->get_data_ptr();
	for (int i = 0; i < in_out->get_number_of_elements(); i++){
		if (x[i] > max) x[i] =  max;

	}
}

/**
 * @brief Calculates the elementwise signum function on the array
 * @param[in,out] in_out Working array
 */
template<class T>
void inplace_sgn(hoNDArray<T> *in_out){
	T* x = in_out->get_data_ptr();
		for (int i = 0; i < in_out->get_number_of_elements(); i++) x[i] = sgn(x[i]);
}

}
