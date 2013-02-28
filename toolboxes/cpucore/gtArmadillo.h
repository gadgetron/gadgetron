#pragma once

#include "hoNDArray.h"
#include <armadillo>

namespace Gadgetron{
template<class T> arma::Mat<typename stdType<T>::type> as_matrix(hoNDArray<T>* x){
	if (x->get_number_of_dimensions() != 2)
		BOOST_THROW_EXCEPTION( runtime_error("Wrong number of dimensions. Cannot convert hoNDArray to matrix"));
	return arma::Mat<typename stdType<T>::type>( (typename stdType<T>::type* ) x->get_data_ptr(),x->get_size(0),x->get_size(1),false);
}


template<class T> arma::Col<typename stdType<T>::type > as_vector(hoNDArray<T>* x){
	return arma::Col< typename stdType<T>::type >((typename stdType<T>::type *) x->get_data_ptr(),x->get_number_of_elements(),false);
}
}
