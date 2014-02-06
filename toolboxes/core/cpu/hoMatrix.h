#pragma once

#include "ho2DArray.h"
#include "complext.h"

#ifdef USE_MKL
    #include "mkl.h"
#endif // USE_MKL

#ifdef GT_Complex8
    #undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
    #undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16; 

namespace Gadgetron{

// the hoMatrix stores every column as the first dimension
// it has the column-wise storage
template <class T> class  hoMatrix : public ho2DArray<T>
{
public:

    typedef hoMatrix<T> Self;
    typedef ho2DArray<T> BaseClass;

    hoMatrix();
    hoMatrix(size_t rows, size_t cols);
    hoMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct = false);

    virtual ~hoMatrix();

    hoMatrix(const hoMatrix<T>& a);
    hoMatrix<T>& operator=(const hoMatrix& rhs);

    virtual bool createMatrix(size_t rows, size_t cols);
    virtual bool createMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t r , size_t c);
    const T& operator()(size_t r , size_t c) const;

    size_t rows() const;
    size_t cols() const;

    // assign the upper/lower triangle matrix as a fixed value
    bool upperTri(const T& v);
    bool lowerTri(const T& v);

    // sum along row or col
    bool sumOverRow(hoNDArray<T>& res) const;
    bool sumOverCol(hoNDArray<T>& res) const;

    // get the sub matrix
    bool subMatrix(Self& res, size_t startR, size_t endR, size_t startC, size_t endC) const;

    bool operator == (const Self& m) const;
    bool operator != (const Self& m) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;
    using BaseClass::accesser_;
};

}

#include <hoMatrix.cpp>
