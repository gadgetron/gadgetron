#pragma once
#ifdef NDEBUG
#define ARMA_NO_DEBUG
#endif


#include "hoNDArray.h"


#include "armadillo"

/** \file hoArmadillo.h
\brief Utilities to create an Armadillo matrix or column vector from an hoNDArray.

Utilities to create an Armadillo matrix or column vector from an hoNDArray.
A helper function that creates an hoNDArray from an Armadillo matrix or vector is deliberatly omitted:
The reccomended approach to using Armadillo's functionality and providing an hoNDArray of the result is 
1) create an hoNDArray to hold the result, 
2) convert this array to an Armadillo matrix or vector using the utilities provided in this header,
3) assign the desired Armadillo computation to this array.
This approach ensures that the Gadgetron -- and not Armadillo -- is responsible for subsequent memory handling.
We refer to hoNDArray_math.h for some specific examples on how to use this Armadillo interface.
*/

namespace Gadgetron{

  /**
   * @brief Creates an Armadillo matrix from a two-dimensional hoNDArray.
   * @param[in] x Input array.
   * @return An Armadillo array mapped to the data pointer of the hoNDArray.
   */
  template<class T> arma::Mat<typename stdType<T>::Type> as_arma_matrix( hoNDArray<T>& x )
  {
    if( x.get_number_of_dimensions() != 2 )
      throw std::runtime_error("Wrong number of dimensions. Cannot convert hoNDArray to matrix");
    return arma::Mat<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_size(0), x.get_size(1), false, true );
  }

  /**
   * @brief Creates an Armadillo matrix from a two-dimensional hoNDArray.
   * @param[in] x Input array.
   * @return An Armadillo array mapped to the data pointer of the hoNDArray.
   */
  template<class T> const arma::Mat<typename stdType<T>::Type> as_arma_matrix( const hoNDArray<T>& x )
  {
    if( x.get_number_of_dimensions() != 2 )
      throw std::runtime_error("Wrong number of dimensions. Cannot convert hoNDArray to matrix");
    return arma::Mat<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_size(0), x.get_size(1), false, true );
  }
  
  /**
   * @brief Creates an Armadillo column vector from an arbitrary-dimensional hoNDArray.
   * @param[in] x Input array.
   * @return An Armadillo array mapped to the data pointer of the hoNDArray.
   */
  template<class T> arma::Col<typename stdType<T>::Type > as_arma_col( hoNDArray<T>& x )
  {
    return arma::Col<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_number_of_elements(), false, true );
  }

  /**
   * @brief Creates an Armadillo column vector from an arbitrary-dimensional hoNDArray.
   * @param[in] x Input array.
   * @return An Armadillo array mapped to the data pointer of the hoNDArray.
   */
  template<class T> const arma::Col<typename stdType<T>::Type > as_arma_col( const hoNDArray<T>& x )
  {
    return arma::Col<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_number_of_elements(), false, true );
  }

  /**
     * @brief Creates an Armadillo row vector from an arbitrary-dimensional hoNDArray.
     * @param[in] x Input array.
     * @return An Armadillo array mapped to the data pointer of the hoNDArray.
     */
    template<class T> arma::Row<typename stdType<T>::Type > as_arma_row( hoNDArray<T>& x )
    {
      return arma::Row<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_number_of_elements(), false, true );
    }

    /**
     * @brief Creates an Armadillo row vector from an arbitrary-dimensional hoNDArray.
     * @param[in] x Input array.
     * @return An Armadillo array mapped to the data pointer of the hoNDArray.
     */
    template<class T> const arma::Row<typename stdType<T>::Type > as_arma_row( const hoNDArray<T>& x )
    {
      return arma::Row<typename stdType<T>::Type>( (typename stdType<T>::Type*) x.get_data_ptr(), x.get_number_of_elements(), false, true );
    }
}


