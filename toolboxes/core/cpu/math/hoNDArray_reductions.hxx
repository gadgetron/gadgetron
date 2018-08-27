#pragma once

#include "cpp_blas.h"
#include "hoNDArray.h"
#include <complex>

namespace Gadgetron {
    template<class T>
     T dot(const hoNDArray<T> *x, const hoNDArray<T> *y, bool cc) {
        return BLAS::dot(x->get_number_of_elements(), x->get_data_ptr(), 1, y->get_data_ptr(), 1);
    }

    template<class T>
     std::complex<T>
    dot(const hoNDArray<std::complex<T>> *x, const hoNDArray<std::complex<T>> *y,
                             bool cc) {
        if (cc) {
            return BLAS::dotc(x->get_number_of_elements(), x->get_data_ptr(), 1, y->get_data_ptr(), 1);
        } else {
            return BLAS::dotu(x->get_number_of_elements(), x->get_data_ptr(), 1, y->get_data_ptr(), 1);
        }
    }


    template<class T>
    T dot (const hoNDArray<T>& x, const hoNDArray<T>& y, bool cc){ return dot(&x,&y,cc);}


    template<class T>
     typename realType<T>::Type asum(const hoNDArray<T> *x) {
        return BLAS::asum(x->get_number_of_elements(), x->get_data_ptr(), 1);
    }

    template<class T>
     typename realType<T>::Type asum(const hoNDArray<T> &x) {
        return asum(&x);
    }


    template<class T>
     size_t amax(const hoNDArray<T> *x) {
        return BLAS::amax(x->get_number_of_elements(), x->get_data_ptr(), 1);
    }

    template<class T>
     size_t amax(const hoNDArray<T> &x) { return amax(&x); }


    template<class T>
     typename realType<T>::Type nrm2(const hoNDArray<T> *x) {
        return BLAS::nrm2(x->get_number_of_elements(), x->get_data_ptr(), 1);
    }

    template<class T>
     typename realType<T>::Type nrm2(const hoNDArray<T>& x) {
        return nrm2(&x);
    }

}


